clc; clear all; close all; 
%% Initial value
load('params.mat');
LINK_NUM = 3;
G = [0.8, 0.03, 0.001; 
    0.03, 0.75, 0.04;
    0.001, 0.04, 0.8];
Gii = diag(G);
P = 0.05*ones(LINK_NUM,1); %inital P
P_new = zeros(LINK_NUM,1);
SINR = zeros(LINK_NUM,1); 
lambda_init = LAGRANIAN_INIT*ones(LINK_NUM,1); %inital lambda
eta1_init = 1; eta1_step = 0.9;
value = [];
%% With Rmin
mu_init = LAGRANIAN_INIT*ones(LINK_NUM,1);
eta2_init = 1; eta2_step = 0.9;
R_min = 0;
%% Gradient descent
mu = mu_init;
lambda = lambda_init; 
eta1 = eta1_init;
eta2 = eta2_init;
for iter_in = 1:MAX_ITER %update lambda
    SINR = Gii.*P./(noise+G*P-Gii.*P);
    val = sum(log2(SINR));
    value = [value;val];
    for i = 1:LINK_NUM %update Pi
        Pi_inner=0;
        for j=1:LINK_NUM
            if j~=i
               Gji = G(j,i);
               Gij = G(i,j);
               Pi_inner = Pi_inner + (1+mu(j))*Gji/(G(:,j)'*P - G(j,j)*P(j)+noise);
            end
        end
        P_new(i) = (1+mu(i))/(lambda(i) + Pi_inner);     
    end
    if abs(P-P_new) < TOLERANCE    
        break;
    end
    eta1 = eta1*eta1_step;
    eta2 = eta2*eta2_step;
    lambda = max(lambda - eta1*(P_max-P_new),0);
    mu = max(mu - eta2*(log2(P_new)-R_min),0);
    P=P_new;
end
P
val
iter_in
figure;
plot(value);