% Author: Hussam Al Daas
% Email: hussam.al-daas@stfc.ac.uk
% Date: 14/04/2022
function Ai = get_local_matrices(A, O)
N = length(O);
Ai = cell(N, 1); % Cell for local problems

for i = 1 : N
    Ai{i} = A(O{i}, O{i});
end