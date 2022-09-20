% Author: Hussam Al Daas
% Email: hussam.al-daas@stfc.ac.uk
% Date: 14/04/2022
function CA = set_up_one_level(A)
% A is a cell that contains the set of local matrices
CA = A; % Cell for Chol of local problems
for i = 1 : length(A)
    try 
        CA{i} = chol(A{i});
    catch
        s = svds(A{i},1);
        CA{i} = chol(A{i} + s*1e-15*speye(length(A{i})));
    end
end