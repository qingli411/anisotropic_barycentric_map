function [h, h2] = plotAnisotropicBarycentricMap2(c, c2, varargin)
% plotAnisotropicBarycentricMap2
%   [h, h2] = plotAnisotropicBarycentricMap2(c,c2) plots the
%   anisotropic barecentric map given the anisotropic barycentric
%   coordinates c and c2. Dot plot(h) for c and line plot(h2) for c2.
%
%   [h, h2] = plotAnisotropicBarycentricMap2(c,c2,wgt,1,cb_label) plots
%   colored map with dots and line colored according to wgt. The last
%   two arguments determine if a colorbar is plotted and the label for
%   the colorbar.
%
%   See also plotAnisotropicBarycentricMap, setupAnisotropicBarycentricMap

    nArgs = length(varargin);
    if nArgs == 0
        l_wgt = 0;
    elseif nArgs == 3
        l_wgt = 1;
        wgt = varargin{1};
        l_colorbar = varargin{2};
        cb_label = varargin{3};
    else
        error(['plotAnisotropicBarycentricMap2(c, c2, '...
        '[Weight, ColorbarOn, CBLabel])'])
    end

    % setup and get the vertices of the barycentric map
    [xc, yc] = setupAnisotropicBarycentricMap();

    % get the Cartesian coordinates
    xx = xc*c';
    yy = yc*c';

    xx2 = xc*c2';
    yy2 = yc*c2';

    % plot data
    if l_wgt
        colormap(jet);
        h = scatter(xx,yy,30,wgt,'o','filled');
        h2 = color_line(xx2,yy2,wgt,'LineWidth',1.5);
        caxis([-1, 0]);
        if l_colorbar
            cb = colorbar('FontSize', 12);
            ylabel(cb,cb_label,'FontSize',14,...
                'Interpreter','tex');
        end
    else
        h = plot(xx,yy,'k.');
        h2 = plot(xx2,yy2,'-k','LineWidth',1.5);
    end

end
