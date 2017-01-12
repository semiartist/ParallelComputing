#include <stdio.h>
#include "mpi.h"
#include "math.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>

/******************************
**DEFINE CALCULATION CONSTANT**
******************************/
const int ROOT = 0;
const double TOL = 0.000000001;
const int ITERATION = 50000;
// const unsigned int N = 10;
const double PI = 3.1415926535897;
/*****************************/
/*******************************************************************************
* Below is the inital calculating process - before enter te loop***************
*******************************************************************************/
void compressed_A_dot_guessX(double **workA, double *workX0, double *workTempB, const unsigned int workStart, const unsigned int workQty, const unsigned int sizeA, double & workNormTempB, const int N){
    /* FUNCTION DOCUMENT
    ** INPUT
    workA = workQty * 5 matrix
    workX0 = sizeA * 1 vector
    workStart = unsigned int, start index
    workQty = unsinged int, row number of workA and workTempB;

    ** OUTPUT
    workTempB = workA * workX0;
    workNormTempB = norm^2 of current node's tempB, need to reduce+ all core;
    */

    // decode a;
    double row;
    workNormTempB = 0;
    unsigned int index = workStart;

    // do the multiplication
    for (size_t i = 0 ; i < workQty ; ++i){
        row = 0;
        int ind0 = index - N;
        int ind1 = index - 1;
        int ind3 = index + 1;
        int ind4 = index + N;

        // test portion
        // printf ("row -> %d, i0 -> %d, i1 -> %d, i3 -> %d, i4 -> %d\n", index, ind0, ind1, ind3, ind4);
        // test portion
        // if (workStart == 15) printf(" --> THE START ROW VALUE IS ->%f\n" , row);

        if (ind0 >= 0 ) row+= workA[i][0] * workX0[ind0];
        if (ind1 >= 0 ) row+= workA[i][1] * workX0[ind1];
        row+= workA[i][2] * workX0[index];
        if (ind3 < sizeA ) row+= workA[i][3] * workX0[ind3];
        if (ind4 < sizeA ) row+= workA[i][4] * workX0[ind4];

        // if (workStart == 15) {
        //     printf(" --> THE ROW VALUE IS ->%f\n" , row);
        //     if (row > 10){
        //         printf("\nROW VALUE TOO BIG:->\n ");
        //     }
        // }
        // WRITE TO workTempB
        workTempB[i] = row;
        // write to workNormTempB
        workNormTempB += row * row;
        // increment index to next row;
        ++index;
    }
}

void get_work_variables(const double *b, const double *workTempB, double *workResidual, double *LRH , double & workNormR, const int workStart, const int workQty, const int N){
    /* FUNCTION DOCUMENT
    ** INPUT
    b = 2N * 1 vector: need to decode first before use;
    workTempB = workQty * 1 vector: current node hold part of tempB,


    ** OUTPUT
    workResidual = workQty * 1 vecor, (part)b - workTempB;
    workH = workResidual/4.  (inv(eye*4)*residual = residual ./ 4)
    workNormR = norm^2 of current node hold residual;
    */
    unsigned int index = workStart; // the index of current R
    unsigned int indexHalf = N * (N-1);
    unsigned int indB = 0;
    workNormR = 0;
    for (size_t i = 0 ; i < workQty ; ++i){
        workResidual[i] = -1 * workTempB[i];
        if (index < N ) {
            indB = index;
            workResidual[i] += b[indB];
        }else if (index >= indexHalf){
            indB = index - indexHalf + N;
            workResidual[i] += b[indB];
        }
        // assembly the output
        workNormR += workResidual[i] * workResidual[i];
        LRH[i+N] = workResidual[i] / 4;

        // workH[i] = workStart + i;
        ++index;
    }  // for loop
    // END OF FUNTION
}

void compressed_A_dot_vector(double **workA,const double * LRH, double *workAh, double &workBeta, const int workStart, const int workQty, const int N){
    /* FUNCTION DOCUMENT
    ** INPUT
    workA = workQty * 5 matrix
    h = sizeA * 1 vector
    workStart = unsigned int, start index
    workQty = unsinged int, row number of workA and workTempB;

    ** OUTPUT
    workAh = workQty * 1;
    workbeta = trans(h)*A*h
    */

    // decode a;
    double row;
    workBeta = 0;
    unsigned int index = N;

    // do the multiplication
    for (size_t i = 0 ; i < workQty ; ++i){
        row = 0;
        int do0 = index - N;
        int do1 = index - 1;
        int do3 = index + 1;
        int do4 = index + N;

        // test portion
        // printf ("row -> %d, i0 -> %d, i1 -> %d, i3 -> %d, i4 -> %d\n", index, ind0, ind1, ind3, ind4);
        // test portion

        row+= workA[i][0] * LRH[do0];
        row+= workA[i][1] * LRH[do1];
        row+= workA[i][2] * LRH[index];
        row+= workA[i][3] * LRH[do3];
        row+= workA[i][4] * LRH[do4];

        // WRITE TO workTempB
        workAh[i] = row;
        // write to workBeta
        workBeta += row * LRH[index];
        // increment index to next row;
        ++index;
    }
}

/******************************************************************************/
using std::endl;
using std::cout;
using std::string;
/********************************MAIN FUNCTION*********************************/
int main(int argc, char * argv[]){

    if (argc == 0){
        cout << " NO INPUT MATRIX SIZE INFORMATION, CHECK YOUR CODE!" << endl;
        return -1;
    }

    unsigned int N = atoi(argv[1]);

    // variables for MPI
    int cpuN, rank, worldSize;
    // variables for calculation
    // pointNumber
    unsigned int sizeB = 2*N , sizeA = N * N, nPoint = N+1;
    MPI_Status status;
    unsigned int workQty, wall;

    /******************************************
    * VARIABLE FOR CALCULATION
    **** GLOBAL VARIABLES
    ******************************************/
    // global matrix <- same through all cores
    double *result, *b;
    int *distWork, *dispH;
    double normB = 0, beta, alpha, gama, lambda, normR = 0, delta;
    int loops = 0;
    bool LOOP = true;


    /******************************************
    * VARIABLE FOR CALCULATION
    **** LOCAL VARIABLES
    ******************************************/
    // local matrix <- different through all core586s
    double *workX0, *workTempB, *workResidual, *workAh, **workA, *LRH;
    double workBeta, workAlpha, workNormTempB, workNormR, workGama;
    unsigned int workStart, workEnd, index= 0 , LRHSize;


    /***************************************************************************
    * Below is the pre processes for the MPI calculating ***********************
    ***************************************************************************/

    // std::cout.setf(ios::showpoint);
    // std::cout.precision(10);
    // std::cout.setf(ios::fixed);

    // initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // ceate a new communicator to do communication
    MPI_Comm myComm;
    int color;
    if(rank < N) color = rank/(int)N;
    else color = MPI_UNDEFINED;

    // printf("my rank -> %d, my color -> %d\n", rank, color);
    // create a new MPI communicator to do the calculation;
    MPI_Comm_split(MPI_COMM_WORLD, color, rank, &myComm);

    // MAKE SURE ONLY VALID CPU RUNS THE PROGRAM
    if (color != MPI_UNDEFINED){

        double tick1, tick2, tick3, tick4, tick5, tick6, reduceTime = 0.0, gathervTime = 0.0, commTime = 0.0; // for the tick from ROOT;
        // string filename = "output" + string(N) + ".txt";
        std::ofstream output("output2k.txt"); // for output file from ROOT;

        // create instance at each core;
        b = new double[sizeB];  // to store the instance of b matrix, as the middle of b matrix are all zeros, so the matrix size is only 2*N
        result = new double[sizeA];  // to store the final result, in ROOT for output.
        // h = new double[sizeA];



        MPI_Comm_size(myComm, &cpuN);  //get the CPU number in the calculating commnucator
        int END = cpuN - 1;

        // create b matrix at ROOT core
        if (rank == ROOT){
            // calculate the  the b vector;
            double h = 1/(double)nPoint;
            double sinValue;
            for(size_t i = 0; i < N ; ++i) {
                sinValue = sin(PI*h*(i+1));
                b[i] = sinValue;
                b[i+N] = sinValue * exp(-1*PI);
                normB += b[i] * b[i];
                normB += b[i+N] * b[i+N];
            }

            // handle the output file;
            output << "The job is for final project of MAE 609 - high performance computing" << std::endl;
            output << "CPU number -> " << cpuN << ", output rank ->" << rank << endl;
            tick1 = MPI_Wtime();

            printf("TOTAL CORE NUMBER -> %d\n", cpuN);

        }

        // get how many rows each core should work with.

        distWork = new int[cpuN];
        dispH = new int[cpuN];
        dispH[0] = 0;
        int workPerCore = sizeA / cpuN;

        // create a warning for the send LRH vector
        // if (workPerCore < N){
        //     cout << "WARNING, THE CODE MAY NOT EXICUTE RIGHT, AS EACH PROCESSOR MAY HANDLE LESS THAN THE N number\n" ;
        // }

        wall = sizeA % cpuN;
        // if(wall == 0) wall += cpuN;
        for (size_t i = 0 ; i<cpuN ; ++i ){
            if (i < wall ) {
                distWork[i] = workPerCore + 1;
            }else {
                distWork[i] = workPerCore;
            }

            // if (i > 0) {
            //     dispH[i] = dispH[i-1] + distWork[i-i];
            //     // printf ("work[i-1] ->%d + disp[i-1]%d = disp[i]%d\n",distWork[i-1], dispH[i-1], dispH[i]);
            //     // printf ("6+5 = %d\n", 6+5);
            // }
            // printf("core -> %d, workload -> %d, dispH -> %d\n", (int)i, distWork[i], dispH[i]);

            // why can not put it inside the for loop?>
            /*
            why
            why
            why
            why
            why
            */
        }
        // printf("core -> 0, workload -> %d, dispH -> %d\n",  distWork[0], dispH[0]);
        for (size_t i = 1; i < cpuN ; ++i) {
            dispH[i] = dispH[i-1] + distWork[i-1];
            // printf("core -> %d, workload -> %d, dispH -> %d\n", (int)i, distWork[i], dispH[i]);
        }

        workQty = distWork[rank];

        // create the LRH matrix;
        LRHSize = workQty + 2*N;
        LRH = new double[LRHSize];

        // make sure every core knows how much they need to work
        // MPI_Scatter(distWork, 1, MPI_INT , &workQty, 1, MPI_INT, 1, MPI_COMM_WORLD);   // send each
        MPI_Bcast(&wall, 1, MPI_INT, 1, myComm);  // broad cast the JUMP of the calculating to all processors
        MPI_Bcast(&normB, 1, MPI_DOUBLE, ROOT, myComm); // broad cast the normB to all processors
        MPI_Bcast(b, 2*N, MPI_DOUBLE, ROOT, myComm); // broad cast b vector;

        // each cpu create their own work matrix;
        // first create workA
        workA = new double*[workQty];  // the chunk A matrix for each processes

        // calculate start and end index of A;
        workStart = workQty * rank;
        if(rank >= wall) workStart += wall;
        index += workStart;

        // generate workA at each core;
        for (size_t i = 0 ; i < workQty ; ++i){
            workA[i] = new double[5];
            if (index >= N) workA[i][0] = -1;
            if (index > 0 && index%N != 0) workA[i][1] = -1;
            workA[i][2] = 4;
            if (((index+1) < sizeA) && ((index+1)%N != 0)) workA[i][3] = -1;
            if (index + N < sizeA) workA[i][4] = -1;
            ++index;
            // if (rank == ROOT) printf("index -> %d\n", index);
        }


        // Then create the guessed X matrix;
        workX0 = new double[workQty];  // initial guess1

        // then create each tempB matrix, residual matrix and h, Ah vectors;
        workTempB = new double[workQty];  // tempB = A * guessX
        workResidual = new double[workQty];  // residual = b - tempB
        // workH = new double[workQty];   // h = inv(I*4) * residual
        workAh = new double[workQty];  // Ah = A * h;

        /***************************************************************************
        * Below is the inital calculating process - before enter the loop***********
        ***************************************************************************/

        // start ticking
        if (rank == ROOT) tick2 = MPI_Wtime();
        // calcualte the A * x0
        // compressed_A_dot_guessX(workA, result, workTempB, workStart, workQty, sizeA, workNormTempB, N);
        workNormTempB = 0;


        // calcualte residual = b - tempb;
        get_work_variables(b, workTempB, workResidual, LRH, workNormR, workStart, workQty, N);

        // calculate the global normR to calculate delta;
        MPI_Allreduce(&workNormR, &normR, 1, MPI_DOUBLE, MPI_SUM, myComm);

        // loop conditions
        delta = normR / normB;
        if (delta < TOL) LOOP = false;

        // MPI_Barrier(MPI_COMM_WORLD);
        while(LOOP && loops < ITERATION){
            // get time before the all reduce;
            if (rank == ROOT) tick5 = MPI_Wtime();

            // gather workH -> h for all processes;
            // MPI_Allgatherv(workH, workQty, MPI_DOUBLE, h, distWork, dispH, MPI_DOUBLE, myComm);
            // alternate for above, just send the left and right's sizeN workH to each other;

            // 1 -> 2
            if (rank % 2 == 0 && rank != END){
                // send to RIGHT
                MPI_Send( &LRH[LRHSize - N-N], N, MPI_DOUBLE,rank+1, 1, myComm);
            }
            if (rank %2 == 1 ){
                //receive from left
                MPI_Recv(&LRH[0], N, MPI_DOUBLE, rank - 1 , 1, myComm, &status);
            }
            // if (loops < 5 ) printf("cpu -> %d FLAG 1 -> REACHED\n", rank);

            // 2->3
            if ( rank % 2 == 1 && rank != END){
                MPI_Send( &LRH[LRHSize - N -N], N, MPI_DOUBLE,rank+1, 2, myComm);
            }
            if (rank %2 == 0 && rank != ROOT){
                MPI_Recv(&LRH[0], N, MPI_DOUBLE, rank - 1 , 2, myComm, &status);
            }
            // if (loops < 5 ) printf("cpu -> %d FLAG 2 -> REACHED\n", rank);

            // 3 ->2
            if ( rank % 2 == 1 && rank != ROOT){
                MPI_Send( &LRH[N], N, MPI_DOUBLE,rank - 1, 3, myComm);
            }
            if (rank %2 == 0 && rank != END){
                MPI_Recv(&LRH[LRHSize - N ], N, MPI_DOUBLE, rank + 1 , 3, myComm, &status);
            }
            // if (loops < 5 ) printf("cpu -> %d FLAG 3 -> REACHED\n", rank);

            //2 -> 1
            if (rank % 2 == 0 && rank != ROOT){
                // send to RIGHT
                MPI_Send( &LRH[N], N, MPI_DOUBLE,rank - 1, 4, myComm);
            }
            if (rank %2 == 1 && rank != END){
                //receive from left
                MPI_Recv(&LRH[LRHSize - N ], N, MPI_DOUBLE, rank + 1 , 4, myComm, &status);
            }
            // if (loops < 5 ) printf("cpu -> %d FLAG 4 -> REACHED\n", rank);


            // get time after the all reduce;
            if (rank == ROOT) {
                tick6 = MPI_Wtime();
                gathervTime += tick6 - tick5;

                // if (loops < 3){
                //
                //     printf("The N -> %d\n",N);
                //     printf("workLoad -> %d\n LRH Size -> %d\n", workQty, LRHSize);
                //
                //     printf("original H vector -> | the LRH h\n");
                //     for (size_t i = 0; i < sizeA ; ++i){
                //         printf("%f   ", h[i]);
                //         if(i < LRHSize) printf (" %f <-%d", LRH[i], i);
                //         printf (" \n");
                //     }
                // }
            }



            // MPI_Barrier(myComm);
            // if (rank == 1 ){
            //
            //     printf("original H vector -> | the LRH h\n");
            //     for (size_t i = 0; i < LRHSize ; ++i){
            //         printf (" %f <-%d \n", workResidual[i], i);
            //     }
            // }

            // MPI_Barrier(myComm);
            // if (rank == 2 && loops < 3){
            //
            //     printf("\nThe CORE -> %d\n",rank);
            //     printf("workLoad -> %d\n LRH Size -> ", workQty, LRHSize);
            //
            //     printf("original H vector -> | the LRH h\n");
            //     for (size_t i = 0; i < sizeA ; ++i){
            //         printf("%f   ", h[i]);
            //         if(i < LRHSize) printf (" %f <-%d", LRH[i], i);
            //         printf (" \n");
            //     }
            // }

            // MPI_Barrier(myComm);
            //
            // if ( loops < 3){
            //
            //     printf("workH -> \n");
            //     for (size_t i = 0; i < workQty ; ++i){
            //         printf (" %f <-%d\n", workH[i], i);
            //         // printf (" \n");
            //     }
            // }

            // MPI_Barrier(myComm);

            // calculate workAh, workBeta
            compressed_A_dot_vector(workA, LRH, workAh, workBeta, workStart, workQty, N);


            // if (loops == 0){
            //     printf("\nThe CORE -> %d\n THE WORKAH MATRIX FOR THIS CORE IS ->\n",rank);
            //     for (size_t i = 0; i < workQty ; ++i){
            //         printf("%f   \n", workAh[i]);
            //     }
            // }




            // if (loops == 0) printf("\ncore -> %d beta -> %f \n", rank, workBeta);
            // calculate workAlpha;
            workAlpha = 0;
            for (size_t i = 0 ; i < workQty ; ++i){
                workAlpha += workResidual[i] * LRH[N + i];
            }


            alpha = 0;
            beta = 0;
            // get time before the all reduce;
            if (rank == ROOT) tick5 = MPI_Wtime();

            // assembly the global alpha
            MPI_Allreduce(&workAlpha, &alpha, 1, MPI_DOUBLE, MPI_SUM, myComm);
            // assembly the global beta;
            MPI_Allreduce(&workBeta, &beta, 1, MPI_DOUBLE, MPI_SUM, myComm);

            // get time after the all reduce;
            if (rank == ROOT) {
                tick6 = MPI_Wtime();
                reduceTime += tick6 - tick5;
            }

            // calculate lambda to update guess_x
            lambda = alpha / beta;
            // if (loops%2 == 0){

            // get time before the all reduce;
            // if (rank == ROOT) tick5 = MPI_Wtime();

            // calcualte the A * x0
            // MPI_Allgatherv(workX0, workQty, MPI_DOUBLE, result, distWork, dispH, MPI_DOUBLE, myComm);
            // get time after the all reduce;
            // if (rank == ROOT) {
            //     tick6 = MPI_Wtime();
            //     gathervTime += tick6 - tick5;
            // }
            // compressed_A_dot_guessX(workA, result, workTempB, workStart, workQty, sizeA, workNormTempB, N);

            // calcualte residual = b - tempb;
            // get_work_variables(b, workTempB, workResidual, LRH , workH, workNormR, workStart, workQty, N);
            workNormR = 0;
            workGama = 0;
            for (size_t i = 0;  i < workQty ; ++i){
                workX0[i] += lambda * LRH[N + i];
                workResidual[i] -= lambda * workAh[i];
                workNormR += workResidual[i] * workResidual[i];
                workGama += workResidual[i]/4 * workAh[i];
            }

            // get time before the all reduce;
            if (rank == ROOT) tick5 = MPI_Wtime();
            MPI_Allreduce(&workGama, &gama, 1, MPI_DOUBLE, MPI_SUM, myComm);

            // calculate the global normR to calculate delta;
            MPI_Allreduce(&workNormR, &normR, 1, MPI_DOUBLE, MPI_SUM, myComm);

            if (rank == ROOT) {
                tick6 = MPI_Wtime();
                reduceTime += tick6 - tick5;
            }

            // if (rank == ROOT ) printf("delta value -> %f, normR value -> %f\n", delta, normR);

            delta = normR / normB;
            if (delta <= TOL) {
                LOOP = false;
                break;
            }

            // update h
            //       h   = inv(eye(4))*    r  -(tran(r)*A*h/beta)* h;
            // | vector |=     I      * vector -     scaler       * v
            // need to calculate (tran(r)*A*h/beta);
            for (size_t i = 0; i < workQty ; ++i) {
                // workH[i] = workResidual[i]/4 + gama/beta * h[workStart + i];
                LRH[i+N] = workResidual[i]/4 + gama/beta * LRH[i+N];
            }
            ++loops;  // increment the loop identicator

        }// major task force of loop;

        // get the result;
        // if( (LOOP == false) && (loops%2 == 0) ){


        if (rank == ROOT) {
            tick3 = MPI_Wtime();
            commTime += gathervTime + reduceTime;
            printf("distWORK   DISPH \n");
            for (size_t i = 0;  i < cpuN ; ++i){
                printf("%d, --  %d\n", distWork[i], dispH[i]);
            }
        }
        printf("workQty -> %d, rank -> %d\n", workQty, rank);
        MPI_Gatherv(workX0, workQty, MPI_DOUBLE, result, distWork, dispH, MPI_DOUBLE, ROOT, myComm);

        // }else {
        //     // workx1 is the result;
        //     MPI_Gatherv( workX0, workQty, MPI_DOUBLE, result, distWork, dispH, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
        // }

        // BELOW ARE OUTPUT PORTION  //

        if (rank == ROOT){
            tick4 = MPI_Wtime();

            output << "FINAL RESULT IS ->" << endl;
            index = 0;
            for (size_t i = 0; i < N ; ++i) {
                for (size_t j = 0; j < N ; ++j){
                    output << " " << result[index] << " ";
                    ++index;
                }
                output << endl;
            }

            cout << endl;

            cout << "CPU AVAILABLE -> " << worldSize<< endl;
            cout << "CPU USED -> " << cpuN << endl;
            cout << "WORKLOAD PER CORE -> " << workQty << " lines" << endl;
            cout << "ITERATION LOOPS ->" << loops << " times"<<endl;
            cout << "FINAL TOLERANCE -> " << delta << endl;
            cout << "TOTAL TIME USED -> " << tick4 - tick1 << endl;
            cout << "LOOP TIME USED -> " << tick3 - tick2 << endl;
            cout << "SEND AND RECEIVE TIME USED -> " << gathervTime << endl;
            cout << "REDUCE ALL TIME USED-> " << reduceTime << endl;
            cout << "COMMUNICATION TIME USED -> " << commTime << endl;

            output.flush();
            output.close();

            // MPI_Barrier(MPI_COMM_WORLD);

            // printf("-- iteration ->%d\n"
            // "-- TOL ->%10f\n"
            // "-- time consume TOTAL -->%f\n"
            // "-- time consume communication -->%10f\n"
            // "-- time comsume calculate -->%f\n"
            // "\n",
            // loops, delta, tick4-tick1, commTime, tick3 - tick2);
        }

        for (size_t i = 0 ; i < workQty ; ++i){
            delete[] workA[i];
        }

        delete[] b;
        delete[] result;
        delete[] workA;
        delete[] workX0;
        delete[] workTempB;
        delete[] workResidual;
        delete[] workAh;
        delete[] distWork;
        delete[] dispH;
        delete[] LRH;
        // if (rank == 1) delete[] distWork;
    } // the if statement for running in myComm

    MPI_Finalize();

    return 0;

}
