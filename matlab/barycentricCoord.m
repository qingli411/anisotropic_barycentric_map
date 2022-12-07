function c = barycentricCoord(a)
% barycentricCoord barycentric coordinates.
%   c = barycentricCoord(a) returns the coordinates of the barycentric map
%   given the anisotropy tensor a_{ij}
%
%   c(1) = lambda(1)-lambda(2);
%   c(2) = 2.*(lambda(2)-lambda(3));
%   c(3) = 3.*lambda(3)+1;
%   where lambda(1:3) is the eigenvalues of a_{ij} in descending order
%
%   See Banerjee et al., 2007

    c = zeros(1,3);
    tmp = eig(a);
    lambda = sort(tmp,'descend');
    c(1) = lambda(1)-lambda(2);
    c(2) = 2.*(lambda(2)-lambda(3));
    c(3) = 3.*lambda(3)+1;    
end % function barycentric_coord