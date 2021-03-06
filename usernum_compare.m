clc; clear all; close all; 
%% Initial value
load('params.mat');
eta1_init = 1; eta1_step = 0.9;
%% With Rmin
eta2_init = 1; eta2_step = 0.9;
R_min = 0;
%% Method
mu = mu_init;
lambda = lambda_init; 
eta1 = eta1_init;
eta2 = eta2_init;
G_num = [5:5:50];
opt_c=zeros(10,2);
opt_r=zeros(10,2);
opt_s=zeros(10,2);

for n = 1:10
    LINK_NUM = G_num(n);
    G = randG(LINK_NUM);
    Gii = diag(G);
    P = unifrnd(0,P_max,LINK_NUM,1); %inital P
    P_new = zeros(LINK_NUM,1);
    SINR = zeros(LINK_NUM,1); 
    lambda_init = LAGRANIAN_INIT*ones(LINK_NUM,1); %inital lambda
    mu_init = LAGRANIAN_INIT*ones(LINK_NUM,1);
    value = [];
    %% With alpha, beta
    SINR = Gii.*P./(noise+G*P-Gii.*P);
    alpha_new = SINR./(1.+SINR);
    beta_new = log(1+SINR)-alpha_new.*log(SINR);
    %% Complicate
    for iter_a_b = 1:MAX_ITER
        alpha = alpha_new;
        beta = beta_new;
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
                       Pi_inner = Pi_inner + (alpha(j)+mu(j))*Gji/(G(:,j)'*P - G(j,j)*P(j)+noise);
                    end
                end
                P_new(i) = (alpha(i)+mu(i))/(lambda(i) + Pi_inner);     
            end
            if abs(P-P_new) < TOLERANCE    
                iter_in
                break;
            end
            eta1 = eta1*eta1_step;
            eta2 = eta2*eta2_step;
            lambda = max(lambda - eta1*(P_max-P_new),0);
            mu = max(mu - eta2*(log(P_new)-R_min),0);
            P=P_new;
        end
        SINR = Gii.*P./(noise+G*P-Gii.*P);
        alpha_new = SINR./(1.+SINR);
        beta_new = log(1+SINR)-alpha_new.*log(SINR);
        if norm(alpha_new-alpha)<TOLERANCE && norm(beta_new-beta)<TOLERANCE
            break
        end
    end
    opt_c(n,1) = val;
    opt_c(n,2) = size(value,1);
    %% With R_min
    value = [];
    P = unifrnd(0,P_max,LINK_NUM,1); %inital P
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
        if abs(P-P_new) < TOLERANCE    
            break;
        end
        eta1 = eta1*eta1_step;
        eta2 = eta2*eta2_step;
        lambda = max(lambda - eta1*(P_max-P_new),0);
        mu = max(mu - eta2*(log(P_new)-R_min),0);
        P=P_new;
    end
    opt_r(n,1) = val;
    opt_r(n,2) = size(value,1);
    %% Simple
    value = [];
    P = unifrnd(0,P_max,LINK_NUM,1); %inital P
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
                   Pi_inner = Pi_inner + Gji/(G(:,j)'*P - G(j,j)*P(j)+noise);
                end
            end
            P_new(i) = 1/(lambda(i) + Pi_inner);     
        end
        if abs(P-P_new) < TOLERANCE

            break;
        end
        eta1 = eta1*eta1_step;
        lambda = max(lambda - eta1*(P_max-P_new),0);
        P=P_new;
    end
    opt_s(n,1) = val;
    opt_s(n,2) = size(value,1);
end
figure;
plot(G_num,opt_c(:,1),'*-r'); hold on
plot(G_num,opt_r(:,1),'o-b'); hold on
plot(G_num,opt_s(:,1),'^-g'); hold on
xlabel('# user'); ylabel('optimum value');
legend('Rmin + alpha','Rmin','Simple','Location','SouthEast');
figure;
plot(G_num,opt_c(:,2),'*-r'); hold on
plot(G_num,opt_r(:,2),'o-b'); hold on
plot(G_num,opt_s(:,2),'^-g'); hold on
legend('Rmin + alpha','Rmin','Simple');
xlabel('# user'); ylabel('iteration times');
