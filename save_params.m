MAX_ITER = 1000;
MAX_ITER_INNER=1000;
LAGRANIAN_INIT = 5000;
TOLERANCE = 1e-5;
P_max = 0.2;
P_min = 0;
G = [0.8, 0.03, 0.001; 
    0.03, 0.75, 0.04;
    0.001, 0.04, 0.8];
save('params.mat')