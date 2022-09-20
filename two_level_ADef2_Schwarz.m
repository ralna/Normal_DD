% Author: Hussam Al Daas
% Email: hussam.al-daas@stfc.ac.uk
% Date: 14/04/2022
function y = two_level_ADef2_Schwarz(A, RAi, loc, RA0, V, x)
Q =@(x) V * (RA0\(RA0'\(V'*x)));
t = Q(x);
t1 = Schwarz(RAi, loc, x);
y = t + t1 - Q(A * t1);