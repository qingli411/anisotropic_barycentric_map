%-------------------------------
% Example
%
% This is an example to visualize the anisotropy of Langmuir turbulence
% as a function of normalized depth using the tools in this repository.
% It reads in and plots the mean Reynolds stress (6 components)
% from a large eddy simulation of Langmuir turbulence
% on an anisotropic barycentric map (figure 1) and the corresponding
% direction of anisotropy map (figure 2). In this example, colorbar shows
% the normalized depth. But if a timeseries of the Reynolds stress is
% given, colorbar may be used to show the time.
%-------------------------------

% close figures and clean up workspace
close all; clear variables;

% load data
data = load('../data/example_langmuir.mat');
uu = data.uu;
vv = data.vv;
ww = data.ww;
uv = data.uv;
uw = data.uw;
vw = data.vw;
z  = data.z;
hb = data.hb;

% find the index at boundary layer base
[~, ind_hb] = min(abs(z+hb));
% weight for color, here fraction in the boundary layer
wgt = -z(1:ind_hb)./z(ind_hb);
% flag for colorbar: 1, add colorbar; 0, no colorbar
l_colorbar = 1;
% colorbar label
cb_label = 'z/h_b';

%-------------------------------
% figure 1: anisotropic barycentric map
%-------------------------------
figure();
setFigProperty;
% compute anisotropy tensor and barycentric coordinate
c = zeros([ind_hb,3]);
for i=1:ind_hb
    a = anisotropyTensor(uu(i), vv(i), ww(i),...
                         uv(i), uw(i), vw(i));
    c(i,:) = barycentricCoord(a);
end
% plot anisotropic barycentric map
plotAnisotropicBarycentricMap(c, wgt, l_colorbar, cb_label);
% save figure
print('anisotropicBarycentricMap', '-dpng', '-r300');

%-------------------------------
% figure 2: direction of anisotropy
%-------------------------------
figure();
setFigProperty;
% compute anisotropy tensor, its eigenvalues and eigenvectors
% corresponding to the maximum and minimum eigenvalue
lambda = zeros([ind_hb,3]);
cmin = zeros([ind_hb,3]);
cmax = zeros([ind_hb,3]);
for i=1:ind_hb
    a = anisotropyTensor(uu(i), vv(i), ww(i),...
                         uv(i), uw(i), vw(i));
    [lambda(i,:),cmax(i,:),cmin(i,:)] = eigMaxMin3(a);
end
% plot direction of anisotropy
plotEigenVectorDirectionMaxMin(lambda, cmax, cmin,...
                               wgt, l_colorbar, cb_label);
% save figure
print('directionOfAnisotropy', '-dpng', '-r300');

