clear all;
close all;
clc;



input500 = [32.829601 , 23.481828, 52.228346, 38.880194, 47.115313, 106.129091, 2837.073885];
vectorSize500 = [36 72 108 144 240 480 1200];

plot(vectorSize500, input500 ,'r-*','LineWidth',2);
xlabel('CPU Number');
ylabel('Time Consumed(s)');
title('matrix size of 500*500');

input1k = [152.546136 193.049571 234.343909 384.852511 1211.192725 ];
vectorSize1000 = [72 144 240 480 1200 ];
figure
plot(vectorSize1000, input1k ,'b-*','LineWidth',2);
xlabel('CPU Number');
ylabel('Time Consumed(s)');
title('matrix size of 1000*1000');

total2k = [162.533 95.0455 73.7426  60.793 61.3212];
comm2k = [ 108.825 74.563  62.8102 53.1454 59.2037];
loop2k = [10869 8879 9670 13964 9164];
cal2k = (total2k - comm2k)./loop2k * 10000;
vectorSize2K = [60 120 240 480 1200]

figure
plot(vectorSize2K, total2k ,'g-*','LineWidth',2);
hold on;
plot(vectorSize2K, cal2k,'r-x','LineWidth',2 )
xlabel('CPU Number');
ylabel('Time Consumed(s)');
title('matrix size of 2000*2000');
legend('total time used','calculation time used');

total5k = [ 1865.72 1233.94 969.878 606.958 354.608 ];
vectorSize5k = vectorSize2K;
figure
plot(vectorSize5k, total5k ,'g-*','LineWidth',2);
xlabel('CPU Number');
ylabel('Time Consumed(s)');
title('matrix size of 5000*5000');
legend('total time used');