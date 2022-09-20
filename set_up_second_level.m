% Author: Hussam Al Daas
% Email: hussam.al-daas@stfc.ac.uk
% Date: 14/04/2022
function [V, A0, lVi] = set_up_second_level(A, Ai, D, Bi, O, tau, k, varargin)

opts.tol = 1e-2;
opts.p = 2 * k;
opts.disp = 0;

for i = 1 : 2 : length(varargin)
    switch varargin{i}
        case 'verbosity'
            opts.disp = varargin{i + 1};
        case 'tolerance'
            opts.tol = varargin{i + 1};
        case 'subspaceDim'
            opts.p = max(k+2, varargin{i + 1});
        case 'saturation'
            saturation =varargin{i + 1};
    end
end
    
n   = length(A);   % Problem size
N   = length(Ai);  % Number of subdomains
lVi = zeros(N, 1); % Numbers of contributed vectors to the 2nd level
W   = cell(N, 1);
for i = 1 : N
    temp = opts;
    DAD = D{i} * (Ai{i} * D{i});
    DAD = (DAD + DAD')/2;
    saturated = 0;
    kk = min(k, length(DAD) - 4);
    try
        rB = chol(Bi{i});
    catch
        s = svds(Bi{i},1);
        rB = chol(Bi{i}+s*1e-14*speye(length(Bi{i})));
    end
    while(saturated == 0)
        temp.p = max(kk+2, min(temp.p, length(DAD) - 2));
        temp.cholB = 1;
        [Z, l] = eigs(DAD, rB, kk, 'largestabs', temp);
        [l,e] = sort(diag(l), 'descend');
        Z = Z(:, e);
        e = find(abs(l) > tau);
        if(length(e) == kk && kk < length(DAD)-4 && saturation == 1)
            kk = min(2*kk, length(DAD)-4);
        else
            saturated = 1;
        end
    end
    Z = D{i} * Z(:, e);
    W{i} = sparse(n, length(e));
    W{i}(O{i}, :) = Z;
    lVi(i) = length(e);
end

slVi = [0; cumsum(lVi)]; % Indices of local contributed vectors to 2nd level
lV   = sum(lVi);         % Size of 2nd level
V    = sparse(n, lV);    % Basis of 2nd level
A0   = sparse(lV, lV);
if lV > 0
    for i = 1 : N
        V(:, slVi(i) + 1 : slVi(i + 1)) = W{i};
    end
    A0 = V' * (A * V);
    A0 = (A0 + A0')/2;
end
