clc; clear all; close all; 
%% Initial value
load('params.mat');
LINK_NUM = 20;
Gii = diag(G);
P = unifrnd(0,P_max,LINK_NUM,1); %inital P
P_new = zeros(LINK_NUM,1);
SINR = zeros(LINK_NUM,1); 
lambda_init = LAGRANIAN_INIT*ones(LINK_NUM,1); %inital lambda
value = [];
%% With Rmin
mu_init = LAGRANIAN_INIT*ones(LINK_NUM,1);
R_min = 0;
%% Gradient descent
mu = mu_init;
lambda = lambda_init; 
eta1 = eta1_init;
eta2 = eta2_init;
for iter_in = 1:MAX_ITER %update lambda
    SINR = Gii.*P./(noise+G*P-Gii.*P);
    val = sum(log(SINR));
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
    eta1 = eta1*eta1_step;
    eta2 = eta2*eta2_step;
    norm(P_max-P_new)
    lambda = max(lambda - eta1*(P_max-P_new),0);
    SINR = Gii.*P_new./(noise+G*P_new-Gii.*P_new);
    mu = max(mu - eta2*(log(SINR)-R_min),0);
    if norm(P-P_new) < TOLERANCE    
        P=P_new;
        break;
    end
    P=P_new;
end
val
figure;
plot(value);
xlabel('Iteration index'); ylabel('Total data rate');
