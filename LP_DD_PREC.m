% Author: Hussam Al Daas
% Email: hussam.al-daas@stfc.ac.uk
% Date: 14/04/2022
function [M1,M2] = LP_DD_PREC(A,T,shift,N,kappa,k,verbosity)

% This function constructs a preconditioner for the linear system arising
% in the IPM method for LP. The linear system matrix is the normal
% equations matrix:  A^T T^-1 A  + shift I
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters:
% A: Constraint matrix. m \times n matrix
% T: Weights. Vector of size m containing strictly positive values
% shift: regularization parameter. Positive scalar
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
% m = 200;
% A = [gallery('poisson',m);
%      spdiags([0*ones(m*m,1) ones(m*m,1) -ones(m*m,1)],-1:1,m*m,m*m)];
% a=10;
% T = 10.^[a-2*a/(2*m*m):-2*a/(2*m*m):-a];
% T=T';
% shift = 1e-8;
% N = 32;
% kappa = 1e2;
% k=100;
% verbosity = 3;
% [M1,M2] = LP_DD_PREC(A,T,shift,N,kappa,k,verbosity);
% b = randn(m*m,1);
% % One level
% [x1,flag1,relres1,iter1,resvec1] = pcg(@(x) A'*(T.\(A*x)), b, 1e-6,100, M1);
% %Two level
% [x2,flag2,relres2,iter2,resvec2] = pcg(@(x) A'*(T.\(A*x)), b, 1e-6,100, M2);
% figure
% semilogy(resvec1,'-b');
% hold on; grid on;
% semilogy(resvec2,'--r');
% legend('ASM','two-level')

partitioning_tic = tic;
[m,n] = size(A);
sT = sqrt(T);
sTm = diag(sparse(1./sT));
A = sTm*A;
if shift > 0
    A = [A;sqrt(shift)*speye(n,n)];
end
%
C = A' * A;
%
[~, Oi, ~, O, ri, D, kc, km] = set_up_DD_normal_sparse(C, A, N, 1);
partitioning_A = toc(partitioning_tic);
%
one_level_tic = tic;
Ci = get_local_matrices(C, O);
RCi = set_up_one_level(Ci);
one_level = toc(one_level_tic);
%
spsd_tic = tic;
Bi = set_up_Neumann_matrices_normal(A, O, ri);
spsd = toc(spsd_tic);
%
if(isempty(kappa) || kappa < 1)
    tau = -1;
    saturation = 0;
else
    tau = (kappa/kc -1)/km;
    saturation = 1;
end
gevp_tic = tic;
[V, A0, lVi] = set_up_second_level(C, Ci, D, Bi, O, tau, k,...
    'verbosity', 0,...
    'subspaceDim', 2*k,...
    'tolerance', 1e-2,...
    'Neumann_Chol', 1, ...
    'saturation', saturation);
gevp = toc(gevp_tic);
%
factor_A0_tic = tic;
try
    RA0 = chol(A0);
catch
    s = svds(A0,1);
    RA0 = chol(A0 + s*1e-15*speye(length(A0)));
end
factor_A0 = toc(factor_A0_tic);
%
if(verbosity > 0)
    fprintf("Partitioning and computing the normal equation: %e\n", partitioning_A);
    fprintf("One level setup: %e\n", one_level);
    fprintf("SPSD splitting computation: %e\n", spsd);
    fprintf("GEVP: %e\n", gevp);
    fprintf("Factoring second level operator: %e\n", factor_A0);
end
if(verbosity > 1)
    for i = 1 : N
        l(i) = length(Oi{i});
    end
    fprintf("Non overlapping subdomain sizes:\n"); disp(l);
    for i = 1 : N
        l(i) = length(O{i});
    end
    fprintf("Overlapping subdomain sizes:\n");  disp(l);
    fprintf("eigenvectors per subdomain:\n");disp(lVi');
    fprintf("dim of coarse space:\n");disp(sum(lVi));
end

M2 =@(x) two_level_ADef2_Schwarz(C, RCi, O, RA0, V, x);
M1 =@(x) Schwarz(RCi,O,x);

