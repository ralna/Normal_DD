clear; clc;
m = 200;
A = [gallery('poisson',m);
    spdiags([0*ones(m*m,1) ones(m*m,1) -ones(m*m,1)],-1:1,m*m,m*m)];
a=10;
T = 10.^[a-2*a/(2*m*m):-2*a/(2*m*m):-a];
T=T';
shift = 1e-8;
N = 32;
kappa = 1e2;
k=100;
verbosity = 3;
[M1,M2] = LP_DD_PREC(A,T,shift,N,kappa,k,verbosity);
b = randn(m*m,1);
% One level
[x1,flag1,relres1,iter1,resvec1] = pcg(@(x) A'*(T.\(A*x)), b, 1e-6,100, M1);
%Two level
[x2,flag2,relres2,iter2,resvec2] = pcg(@(x) A'*(T.\(A*x)), b, 1e-6,100, M2);
figure
semilogy(resvec1,'-b');
hold on; grid on;
semilogy(resvec2,'--r');
legend('ASM','two-level')