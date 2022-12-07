function [lambda, cmax, cmin] = eigMaxMin3(a)
% eigMaxMin3
%   [lambda, cmax, cmin] = eigMaxMin3(a) returns the eigenvalues lambda
%   and eigenvectors corresponding to the maximum and minimum eigenvalue
%   of a 3x3 matrix
%
%   See also eig

    [E,D] = eig(a);     
    F = [D(1,1) D(2,2) D(3,3);E];
    FF = sortrows(F',-1);
    cmax = FF(1,2:4)';
    cmin = FF(3,2:4)';
    lambda = FF(1:3,1);
end