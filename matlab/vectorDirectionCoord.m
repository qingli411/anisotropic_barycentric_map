function [x, y, s] = vectorDirectionCoord(c, r)
% vectorDirectionCoord
%   [x, y, s] = vectorDirectionCoord(c, r) returns the x- and
%   y-coordinates of a unit vector c in the vector direction map.
%   r is the length of vector, returned from setupVectorDirection().
%   s=1 if c is directed in the positive z-direction, s=-1 otherwise.
%
%   See also setupVectorDirection

    c = c.*r;
    rc = sqrt(c(1)^2+c(2)^2);
    phi = atan2(c(3), rc);
    theta = atan2(c(2), c(1));
    rr = 1-abs(phi/pi*2);
    x = rr*cos(theta);
    y = rr*sin(theta);
    s = sign(phi);
end