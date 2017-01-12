% THIS SCRIPT IS TO PLOT THE RESULT FOR THE FINAL PROJECT OF
% MAE 609 - HIGH PERFORMANCE COMPUTING
% FEI CHEN
% fchen29@buffalo.edu
% 12/20/2016

clear all;
close all;
clc;

n = 100;
increment = 1/n;

% realA matrix is the real value of the matrix
realA = zeros(n+1, n+1);

for i = 1:n+1;
    for j = 1:n+1;
        realA(i,j) = sin(pi*increment*(j-1))*exp(-pi* (increment*(i-1)));
    end
end

%matrix A contains the 200*200 matrix of result;
load ('a.mat');

h = 0:increment:1;

figure(1);
h2 = 0:(1/(size(A,1) + 1)):1;
plot(h(2:end-1) , realA(end-1,2:end-1),'r-x', 'LineWidth' ,2);
hold on;
plot(h2(2:end-1) , A(end,:), 'b-.', 'LineWidth' ,2 );
axis equal;
xlim([0, 1]);
ylim([-0.1, 0.2]);
legend('real Value','calculation result');

figure(2);
plot(h(2:end-1) , realA(50,2:end-1),'r-x', 'LineWidth' ,2);
hold on;
plot(h2(2:end-1) ,  A(99,:), 'b-.', 'LineWidth' ,2 );
axis equal;
legend('real Value','calculation result');


figure(3);
plot(h(2:end-1) , realA(2,2:end-1),'r-x', 'LineWidth' ,2);
hold on;
plot(h2(2:end-1) ,  A(2,:), 'b-.', 'LineWidth' ,2 );
axis equal;
legend('real Value','calculation result');

figure(5);
plot(h(2:end-1) , realA(2:end-1,51),'r-x', 'LineWidth' ,2);
hold on;
plot(h2(2:end-1) ,  A(:,101), 'b-.', 'LineWidth' ,2 );
axis equal;
legend('real Value','calculation result');

figure(4)
surf(h(2:end-1), h(2:end-1), realA(2:end-1,2:end-1));
surf(h2(2:end-1), h2(2:end-1), A);
% legend('real Value','calculation result')