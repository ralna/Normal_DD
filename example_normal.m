clear; clc;
m = 200;
A = [gallery('poisson',m);
     spdiags([0*ones(m*m,1) ones(m*m,1) -ones(m*m,1)],-1:1,m*m,m*m)];
N = 32;
kappa = 5e1;
k=100;
verbosity = 3;
[M1,M2] = NORMAL_DD_PREC(A,N,kappa,k,verbosity);
b = randn(m*m,1);
% One level
[x1,flag1,relres1,iter1,resvec1] = pcg(@(x) A'*(A*x), b, 1e-6,100, M1);
%Two level
[x2,flag2,relres2,iter2,resvec2] = pcg(@(x) A'*(A*x), b, 1e-6,100, M2);
figure
semilogy(resvec1,'-b');
hold on; grid on;
semilogy(resvec2,'--r');
legend('ASM','two-level')