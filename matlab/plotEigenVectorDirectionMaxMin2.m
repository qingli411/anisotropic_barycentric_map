function h = plotEigenVectorDirectionMaxMin2(lambda, cmax, cmin,...
               lambda2, cmax2, cmin2, wgt, l_colorbar, cb_label)
% plotEigenVectorDirectionMaxMin2
%   h = plotEigenVectorDirection2(lambda, cmax, cmin,...
%       lambda2, cmax2, cmin2, wgt, l_colorbar, cb_label)
%   plots the directions of the dominate direction of the turbulence
%   anisotropy.
%
%   The turbulence anisotropy is described by the eigenvalue lambda and
%   the eigenvectors associated with the greatest and smallest eigenvalue
%   cmax and cmin of the anisotropy tensor a_{ij}. If the turbulence
%   anisotropy is cigar-like (near one component), the direction of
%   cmax is marked by O; if the turbulence anisotropy is pancake-like
%   (near two component), the direction of cmin is marked by X.
%
%   Two sets of data are plotted. Set 1 is colored according to wgt.
%   Set 2 is plotted in gray. The flag l_colorbar determines
%   if a colorbar is plotted. String cb_label is the colorbar label.
%
%   Note that the direction of eigenvectors can have both signs. This map
%   shows the direction with positive z.
%
%   See also setupVectorDirection, vectorDirectionCoord,
%   plotEigenVectorDirectionMaxMin

    % set up the vector direction plot
    r = setupVectorDirection();

    % find the coordinates
    [tmpxmax, tmpxmin, tmpymax, tmpymin] = findCoord(cmax, cmin, r);
    [tmpxmax2, tmpxmin2, tmpymax2, tmpymin2] = findCoord(cmax2, cmin2, r);

    % plot data
    colormap(jet);
    bsize = 200;
    % gray for group 2
    cgray = [0.5, 0.5, 0.5];
    eta1 = lambda2(:,1)-lambda2(:,2);
    eta2 = 2.*(lambda2(:,2)-lambda2(:,3));
    inds1 = find(eta1>=eta2);
    inds2 = find(eta1<eta2);
    scatter(tmpxmax2(inds1),tmpymax2(inds1),...
        eta1(inds1).*bsize,'o',...
        'MarkerEdgeColor',cgray,...
        'LineWidth',1.5);
    scatter(tmpxmin2(inds2),tmpymin2(inds2),...
        eta2(inds2).*bsize,'x',...
        'MarkerEdgeColor',cgray,...
        'LineWidth',1.5);
    % color for group 1
    eta1 = lambda(:,1)-lambda(:,2);
    eta2 = 2.*(lambda(:,2)-lambda(:,3));
    inds1 = find(eta1>=eta2);
    inds2 = find(eta1<eta2);
    scatter(tmpxmax(inds1),tmpymax(inds1),...
        eta1(inds1).*bsize,wgt(inds1),'o',...
        'LineWidth',1.5);
    scatter(tmpxmin(inds2),tmpymin(inds2),...
        eta2(inds2).*bsize,wgt(inds2),'x',...
        'LineWidth',1.5);
    scatter(r,r.*1.1,bsize,'o','LineWidth',1.5,...
        'MarkerEdgeColor','k');
    scatter(r,r.*0.9,bsize,'x','LineWidth',1.5,...
        'MarkerEdgeColor','k');
    caxis([-1,0]);
    if l_colorbar
        cb = colorbar('FontSize',12);
        ylabel(cb,cb_label,'FontSize',14,...
                'Interpreter','tex');
    end
    h = gcf;
end

function [tmpxmax, tmpxmin, tmpymax, tmpymin] = findCoord(cmax, cmin, r)
    % check size of c
    csize = size(cmin);
    nc = csize(1);
    xxmin = zeros(nc,1);
    yymin = zeros(nc,1);
    zzmin = zeros(nc,1);
    xxmax = zeros(nc,1);
    yymax = zeros(nc,1);
    zzmax = zeros(nc,1);

    % get the coordinates
    for i=1:nc
        [xxmin(i), yymin(i), zzmin(i)] = ...
            vectorDirectionCoord(cmin(i,:), r);
        [xxmax(i), yymax(i), zzmax(i)] = ...
            vectorDirectionCoord(cmax(i,:), r);
    end
    % adjust the coordinates
    [tmpxmin, tmpymin] = adjustXY(xxmin, yymin, zzmin);
    [tmpxmax, tmpymax] = adjustXY(xxmax, yymax, zzmax);
end

function [tmpx, tmpy] = adjustXY(xx, yy, zz)
    % adjust the coordinates to show the direction with positive z
    tmpx = sign(zz).*xx;
    tmpy = sign(zz).*yy;
end
