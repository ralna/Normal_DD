% Author: Hussam Al Daas
% Email: hussam.al-daas@stfc.ac.uk
% Date: 14/04/2022
function y = Schwarz(RAi, loc, x)
N = size(RAi, 1);
y = 0*x;
for i = 1 : N
    y(loc{i}, :) = y(loc{i}, :) + (RAi{i}\(RAi{i}'\x(loc{i}, :)));
end