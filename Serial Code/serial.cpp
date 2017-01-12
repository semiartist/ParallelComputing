#include <iostream>
#include <stdio.h>
// #include "mpi.h"
#include "math.h"

void generate_A_and_b(unsigned int, double **, unsigned int , double ** , double *);
void cgSolver(double **, unsigned int, double **, double *, unsigned int);
void matrix_dot_vector(double **, double *, double *, unsigned int);

#define N 5
#define PI 3.141592654
#define ITERATION 1000
#define OUTPUTSIZE 5

int main (int argc, char *argv[]){
    // code of serial part;
    // assume the X direction and Y direction has a same finit node;

    unsigned int BSize = (N-1)*(N-1);
    unsigned int ASize = N+1;

    double ** A = new double*[ASize];
    double ** B = new double*[BSize];
    for (size_t i = 0; i < BSize ; ++i){
        B[i] = new double[BSize];
        if (i < ASize) A[i] = new double[ASize];
    }
    double *b = new double[BSize];

    // first assembly the large matrix with boundry condition;
    double h = 1/(double)N;
    for (size_t i = 0; i < ASize ; ++i){
        double sinValue = sin(PI*h*i);
        A[0][i] = sinValue;
        A[ASize-1][i] = sinValue * exp(-1*PI);
    }

    // Then assembly the matrix for solving the problem;
    generate_A_and_b(ASize, A, BSize, B, b);
    // test output the A B b matrix;
    printf("A Matrix is: \n");
    for (size_t i = 0; i < OUTPUTSIZE ; ++i){
        for (size_t j = 0; j < OUTPUTSIZE ; ++j) {
            printf("%f  ",A[i][j]);
        }
        printf("<- line %d\n",(int)i);
    }

    printf("\nB Matrix is: \n");
    for (size_t i = 0; i < OUTPUTSIZE ; ++i){
        for (size_t j = 0; j < OUTPUTSIZE ; ++j) {
            printf("%f  ",B[i][j]);
        }
        printf("<- line %d\n",(int)i);
    }
    printf("\nb vector is: \n");
    for (size_t j = 0; j < OUTPUTSIZE ; ++j) {
        printf("%f\n",b[j]);
    }

    // solve the Ax = b Problem with CG method with Jacobi Preconditioner;
    cgSolver(B, BSize, A, b, ASize);

    printf("After calculation : A Matrix is: \n");
    for (size_t i = 0; i < OUTPUTSIZE ; ++i){
        for (size_t j = 0; j < OUTPUTSIZE ; ++j) {
            printf("%f  ",A[i][j]);
        }
        printf("<- line %d\n",(int)i);
    }



    // clean up;
    for (size_t i = 0; i < BSize ; ++i){
        delete[] B[i];
        if (i < ASize) delete[] A[i];
    }
    delete[] A, B, b;

    /*
    // code of paralle part; TODO
    int i, myid, nprocs, namelen;
    char name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);

    /* get the processor numbers; TODO */

    /* use the rank 0 one to produce the matrix, and send it to
    * the other processors TODO */


    /* Gather information from all processors TODO */


    /* Output result to a file for further analysis TODO */

    // MPI_Finalize();


    return (0);
}

void generate_A_and_b (unsigned int sizeA, double **A, unsigned int sizeB, double **B, double *b ){
    for (size_t i = 1; i< sizeA-1 ; ++i){
        b[i-1] = A[0][i];
        b[sizeB-sizeA+i+1] = A[sizeA-1][i];
    }

    for (size_t i = 0; i < sizeB ; ++i){
        B[i][i] = 4;
        if ((i>0) && ((i-1)%sizeB !=0)) B[i][i-1] = -1;
        if (i>N-2) B[i][i-N+1] = -1;
        if ((i+1 < sizeB) &&(i%sizeB!=0) ) B[i][i+1] = -1;
        if (i+N-1 < sizeB ) B[i][i+N-1] = -1;

    }

}

void cgSolver(double ** B, unsigned int sizeB, double ** A, double * b, unsigned int sizeA){
    // create a Jacobi;

    /*********************
    ****  Initialzing ****
    **********************/

    // DO THE PRECONDITIONING, SIMPLY TIME THE WHOLE DIAG BY 4, AND THE b BY 4
    // create initial guess and set value to 0;
    double * iniX = new double[sizeB];
    double * nextX = new double[sizeB];
    // create other matrix for calculating;
    double * residual = new double[sizeB];
    double * h = new double[sizeB];
    double * newb = new double[sizeB];
    double * a = new double[sizeB];
    double * tau = new double[sizeB];
    // create variables for calculating;
    double alpha, beta, delta, lambda, normb;
    // step 1
    for (size_t i = 0; i < sizeB ; ++i) {
        iniX[i] = 0;
        nextX[i] = 0;
    }


    /*********************
    ** Done Initialzing **
    **********************/

    // step 2
    matrix_dot_vector(B, iniX, newb, sizeB);
    // step 3 and 4
    for (size_t i = 0 ; i < sizeB ; ++i){
        residual[i] =  b[i] - newb[i];
        h[i] = residual[i]*(0.25);
    }
    alpha = 0;
    normb = 0;
    for (size_t i = 0; i < sizeB ; ++i){
        alpha += residual[i]*residual[i];
        normb += b[i] * b[i];
    }
    // step 5
    delta = alpha / normb;

    int it = 0;
    while (delta > 0.00001 && it < ITERATION){
        // step 7
        matrix_dot_vector(B, h, a, sizeB);
        // printf("\n -> computation flag <- \n");
        beta = 0;
        alpha = 0;
        // step 8 & 9
        for (size_t i = 0; i < sizeB ; ++i){
            beta += a[i] * h[i];
            alpha += residual[i]*h[i];  // residual * h / buffer * h == trans(h)*B*h
        }
        // step 10
        lambda = alpha / beta;
        // step 11
        if(it%2 ==0){  // first iteration and so in this loop will store x in next
            for (size_t i = 0; i < sizeB ; ++i){
                nextX[i] = iniX[i] + lambda*h[i];
            }
            matrix_dot_vector(B, nextX, newb, sizeB);
        } else {  // second iteration and so will sotre x in the iniX;
            for (size_t i = 0; i < sizeB ; ++i){
                iniX[i] = nextX[i] + lambda*h[i];
            }
            matrix_dot_vector(B, iniX, newb, sizeB);
        }
        alpha = 0;
        // step 12 & 12+
        for (size_t i = 0; i < sizeB ; ++i){
            residual[i] = residual[i] - lambda * a[i];
            alpha += residual[i] * residual[i];
            // step 13
            tau[i]=  residual[i] / 4;
            // direction[i] = -1 * h[i];
        }
        // step 12+
        delta = alpha / normb;
        if(delta < 0.00001) break;

        alpha = 0;
        for (size_t i = 0 ; i < sizeB ; ++i){
            alpha += tau[i] * a[i];
        }
        double gama = alpha / beta;
        for (size_t i = 0; i < sizeB ; ++i){
            h[i] = tau[i] - gama*h[i];
        }

        printf("\nh vector is: \n");
        for (size_t j = 0; j < sizeB ; ++j) {
            printf("%f\n",h[j]);
        }

        ++it; // inrement the iteation;
    }

    std::cout << "it ->" << it << std::endl;
    std::cout << "delta ->" << delta << std::endl;
    unsigned int ind = 0;
    if (it%2 == 0){
        for (size_t i = 1; i < (sizeA - 1); ++i){
            for (size_t j = 1; j < (sizeA-1); ++j){
                A[i][j] = iniX[ind++];
            }
        }
    } else {
        for (size_t i = 1; i < (sizeA - 1); ++i){
            for (size_t j = 1; j < (sizeA-1); ++j){
                A[i][j] = nextX[ind++];
            }
        }
    }


    // SOLVE THE PROBLEM WITH CONJUGATE DIRECTION METHOD;

    delete[] residual, h, a, tau, iniX, nextX, newb;
}

void matrix_dot_vector(double ** A, double * b, double * c, unsigned int size){
    double temp;
    for (size_t i = 0; i < size ; ++i){
        temp = 0.0;
        for(size_t j = 0; j < size ; ++j){
            temp += A[i][j] * b[j];
        }
        c[i] = temp;
    }

    return;
}

double getNorm(unsigned int vSize, double vector[]){
    double result;
    for (size_t i = 0; i < vSize ; ++i){
        result += vector[i]*vector[i];
    }

    return result;
}
