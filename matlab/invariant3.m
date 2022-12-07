function [I, II, III] = invariant3(a)
% invariant3 invariants of a 3x3 tensor
%   [I, II, III] = invariant3(a) returns the three invariants
%   of a 3x3 tensor:
%
%   I   = a_{ii}
%   II  = a_{ij} a_{ji}
%   III = a_{ij} a_{jk} a_{ki}

%     I = trace(a);
    I = a(1,1)+a(2,2)+a(3,3);
    II = 0;
    III = 0;
    for i=1:3
        for j=1:3
            II = II + a(i,j).*a(j,i);
            for k=1:3
                III = III + a(i,j).*a(j,k).*a(k,i);
            end
        end
    end
end