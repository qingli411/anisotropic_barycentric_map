function h = plotAnisotropicBarycentricMap(c, varargin)
% plotAnisotropicBarycentricMap
%   h = plotAnisotropicBarycentricMap(c) plots the
%   anisotropic barecentric map given the anisotropic barycentric
%   coordinates c.
%
%   h = plotAnisotropicBarycentricMap(c,wgt,1,cb_label) plots colored
%   map with dots colored according to wgt. The last two arguments
%   determine if a colorbar is plotted and the label for colorbar.
%
%   See also plotAnisotropicBarycentricMap2, setupAnisotropicBarycentricMap

    nArgs = length(varargin);
    if nArgs == 0
        l_wgt = 0;
    elseif nArgs == 3
        l_wgt = 1;
        wgt = varargin{1};
        l_colorbar = varargin{2};
        cb_label = varargin{3};
    else
        error(['plotAnisotropicBarycentricMap(c, ',...
            '[Weight, ColorbarOn, CBLabel])']);
    end

    % setup and get the vertices of the barycentric map
    [xc, yc] = setupAnisotropicBarycentricMap();

    % get the Cartesian coordinates
    xx = xc*c';
    yy = yc*c';

    % plot data
    if l_wgt
        colormap(jet);
        h = scatter(xx,yy,30,wgt,'o','filled');
        caxis([-1, 0]);
        if l_colorbar
            cb = colorbar('FontSize', 12);
            ylabel(cb,cb_label,'FontSize',14,...
                'Interpreter','tex');
        end
    else
        h = plot(xx,yy,'k.');
    end

end
