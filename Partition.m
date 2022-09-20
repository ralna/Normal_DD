% Author: Hussam Al Daas
% Email: hussam.al-daas@stfc.ac.uk
% Date: 14/04/2022
function [p1, edgecut, s, beginIn, endIn] = Partition(A, part)
%DefineAlpha performs kway partitioning using metis and defines the parts for each core

%Input:
% A: a square sparse matrix
% part: number of partitions

%Output:
% p1: the permutation vector of metismex(kway)
% edgecut: the number of edges that are split between different domains
%       (output of metismex)
% s: the sizes of each of the subdomains
% beginIn: the starting indices of the subdomains
% endIn: the ending indices of the subdomains
M = A + A';
M = M - diag(diag(M));
cc = length(A);
if (cc == part)
    for i = 1 : part
        beginIn(i) = i;
    end
    endIn = beginIn;
    p1 = beginIn;
    edgecut = 0;
    s = 0;
else
    [parti,edgecut] = metismex('PartGraphKway',M,part, 0);
    xx = 1;
    for i = 1 : part
        bb = find(parti == i - 1);
        p1(xx : xx + length(bb) - 1) = bb;
        beginIn(i) = xx;
        s(i) = length(bb);
        xx = xx + length(bb);
        endIn(i)= xx - 1;
    end  
end
end

