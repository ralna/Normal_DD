% Author: Hussam Al Daas
% Email: hussam.al-daas@stfc.ac.uk
% Date: 14/04/2022
function [p, Oi, Og, O, ri, D, kc, km] = set_up_DD_normal_sparse(C, A, N, zeroBoundaryPartUnity, varargin)
[m, n] = size(C);
if m ~= n
    fprintf("Matrix is not squared\n");
    return;
end

[p, ~, ~, b, e] = Partition(C, N);
Oi = cell(N, 1); % Interior subdomain
Og = Oi;         % Interface of subdomain
O  = Og;         % Overlapping subdomain
D  = O;          % Partition of unity (PoU)
ri = D;          % Nonzero rows of interior subdomain columns in A

p1n = 1:n;
p = p1n(p);
for i = 1 : N
    Oi{i} = p(b(i) : e(i));
    [~, ci, ~]    = find(C(Oi{i}, :));   % Find nonzero column in row block
    ci            = unique(ci)';         % Get the unique ones: this is the overlapping subdomain not ordered
    Og{i}         = setdiff(ci, Oi{i});  % Set the boundary
    O{i}          = [Oi{i}, Og{i}];      % Set the overlapping nodes in order [Interior, Interface]
    [ri{i}, ~, ~] = find(A(:, Oi{i}));   % Get indices of nonzero rows in interior subdomain columns in A
    ri{i} = unique(ri{i});
    D{i}          = speye(length(O{i})); % Set the PoU matrix to I
end

SD = sparse(n, n); % will contain the sum of the partition of unity
for i = 1 : N
    GD             = sparse(n, n); % To expand the local PoU
    GD(O{i}, O{i}) = D{i};         % Set the values
    
    SD = SD + GD; % Sum the ith PoU
end

% Set up linear continuous PoU
for i = 1 : N
    D{i} = sparse(diag(1./diag(SD(O{i}, O{i}))));
end

km = max(diag(SD)); % multiplicity constant

if(zeroBoundaryPartUnity)
    for i = 1 : N
        lOi  = length(Oi{i});
        lO   = length(O{i});
        D{i} = speye(lO);
        D{i}(lOi + 1 : lO, lOi + 1 : lO) = 0 * D{i}(lOi + 1 : lO, lOi + 1 : lO);
    end
end

k_c = zeros(N, 1); % number of colors

for i = 1 : N
    for j = 1 : N
        if j == i
            continue;
        end
        if(~isempty(intersect(O{i}, O{j})))
            k_c(i) = k_c(i) + 1;
        end
    end
end
kc = max(k_c);
kcp1 = kc + 1;
kc = kcp1;