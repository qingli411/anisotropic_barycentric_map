 function a = anisotropyTensor(u1u1,u2u2,u3u3,u1u2,u1u3,u2u3)
% anisotropyTensor turbulent anisotropy tensor.
%   a = anisotropyTensor(u1u1,u2u2,u3u3,u1u2,u1u3,u2u3) returns a 3x3 
%   matrix a_{ij} measuring the anisotropy of the turbulence.
%
%   a_{ij}=\frac{\overline{u_i u_j}}{2TKE}-\delta_{ij}/3,
%   TKE=\frac{\overline{u_i u_i}}{2}
%
%   See more in Lumley, 1979

    a = zeros(3,3);
    tke = 0.5.*(u1u1+u2u2+u3u3);
    a(1,1) = 0.5.*u1u1./tke-1./3;
    a(2,2) = 0.5.*u2u2./tke-1./3;
    a(3,3) = 0.5.*u3u3./tke-1./3;
    a(1,2) = 0.5.*u1u2./tke;
    a(1,3) = 0.5.*u1u3./tke;
    a(2,1) = 0.5.*u1u2./tke;
    a(2,3) = 0.5.*u2u3./tke;
    a(3,1) = 0.5.*u1u3./tke;
    a(3,2) = 0.5.*u2u3./tke;
end % function anisotropy_tensor

