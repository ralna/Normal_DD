% Author: Hussam Al Daas
% Email: hussam.al-daas@stfc.ac.uk
% Date: 14/04/2022
function [M1,M2] = NORMAL_DD_PREC(A,N,kappa,k,verbosity)

% This function constructs a preconditioner for the normal equation linear 
% system. The linear system is:  A^T A x = b
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters:
% A: sparse matrix. m \times n matrix
% N: Number of subdomain. Integer >= 1
% kappa: Upper bound on the condition number. Must be > 1.
% k: Number of eigenvectors to be selected in each subdomain to be 
% included in the second level. Positive integer.
% verbosity: level of printing information.
%           =0: No printing
%           >=1: Print runtime for different stages of the construction
%           >=2: Print subdomain sizes, number of eigenvalues selected per 
%                subdomain
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output parameters:
% M1: a function handle that applies the one-level preconditioner to a 
% vector.
% M2: a function handle that applies the two-level preconditioner to a 
% vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Notes:
% The user can provide kappa and/or k. If both are provided, k is ignored.
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example:
% clear; clc;
% m = 200;
% A = [gallery('poisson',m);
%      spdiags([0*ones(m*m,1) ones(m*m,1) -ones(m*m,1)],-1:1,m*m,m*m)];
% N = 32;
% kappa = 5e1;
% k=100;
% verbosity = 3;
% [M1,M2] = NORMAL_DD_PREC(A,N,kappa,k,verbosity);
% b = randn(m*m,1);
% % One level
% [x1,flag1,relres1,iter1,resvec1] = pcg(@(x) A'*(A*x), b, 1e-6,100, M1);
% %Two level
% [x2,flag2,relres2,iter2,resvec2] = pcg(@(x) A'*(A*x), b, 1e-6,100, M2);
% figure
% semilogy(resvec1,'-b');
% hold on; grid on;
% semilogy(resvec2,'--r');
% legend('ASM','two-level')

shift = 0;
m = size(A,1);
T = ones(m,1);
[M1,M2] = LP_DD_PREC(A, T, shift, N, kappa, k, verbosity);