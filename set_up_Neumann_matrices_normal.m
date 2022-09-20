% Author: Hussam Al Daas
% Email: hussam.al-daas@stfc.ac.uk
% Date: 14/04/2022
function Bi = set_up_Neumann_matrices_normal(A, O, ri)
N = length(O);
Bi = cell(N, 1);
for i = 1 : N
    normai = (svds(A(ri{i}, O{i}),1))^2;
    if(isnan(normai))
        normai = 1;
    end
    loc_shift = 1e-15*normai;
    Bi{i} = A(ri{i}, O{i})' * A(ri{i}, O{i}) + loc_shift * speye(length(O{i}));
end
