MAX_ITER = 1000;
MAX_ITER_INNER=10000;
LAGRANIAN_INIT = 1;
TOLERANCE = 1e-5;
P_max = 0.2;
P_min = 0;
noise = 1e-15;
eta1_init = 1; eta1_step = 1;
eta2_init = 1; eta2_step = 1;
load('G.mat'); 
save('params.mat')