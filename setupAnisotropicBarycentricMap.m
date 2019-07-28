function [xc, yc]  = setupAnisotropicBarycentricMap(varargin)
% setupAnisotropicBarycentricMap
%   [xc, yc] = setupAnisotropicBarycentricMap() sets up the anisotropic
%   barycentric map and returns the x- and y-coordinates of the three
%   vertices: (0.5, 0), (-0.5, 0) and (0, sqrt(3)/2). No arguments are
%   required.
%
%   xc = [0.5, -0.5,           0];
%   yc = [   0,   0, sqrt(3)*0.5];
%
%   So the x- and y-coordinates xx(1,N) and yy(1,N) are related with
%   the barycentric coordinates c(N,3) by
%   xx = xc*c';
%   yy = yc*c';
%
%   See also barycentricCoord

    nArgs = length(varargin);
    if nArgs > 0
        error('No arguments required.');
    end

    % Vertices of the barycentric map
    xc = [0.5, -0.5,           0];
    yc = [   0,   0, sqrt(3)*0.5];

    % plot the triangle
    plot(xc, yc, 'k.');
    hold on;
    for i=1:3
        ip1 = mod(i,3)+1;
        line([xc(i),xc(ip1)],[yc(i),yc(ip1)],...
             'Color', 'black',...
             'LineWidth', 2);
    end

    % grid
    nsp = 5;
    lc = abs(xc(2)-xc(1));
    dc = lc/nsp;
    cl = zeros(3*(nsp-1),3);
    for i=1:3
        ip1=mod(i,3)+1;
        for j=1:nsp-1
            k = (i-1)*(nsp-1)+j;
            cl(k,i) = dc*j;
            cl(k,ip1) = dc*(nsp-j);
        end
    end
    xl = xc*cl';
    yl = yc*cl';

    nl = length(cl);
    for i=1:3
        for j=1:nsp-1
            k = (i-1)*(nsp-1)+j;
            kp = mod(i*(nsp-1)+nsp-j-1,nl)+1;
            line([xl(k),xl(kp)],[yl(k),yl(kp)],...
                 'LineStyle', '--',...
                 'Color', 'black',...
                 'LineWidth', 1);
        end
    end

    % plane strain limit
    c1_ps = [1/3,2/3,0];
    x1_ps = xc*c1_ps';
    y1_ps = yc*c1_ps';
    c2_ps = [0,0,1];
    x2_ps = xc*c2_ps';
    y2_ps = yc*c2_ps';
    line([x1_ps,x2_ps],[y1_ps,y2_ps],...
         'Color', 'black',...
         'LineWidth', 1);

    % turn off axis
    axis off;
    % set aspect ratio to 1:1:1
    daspect([1,1,1]);

    % labels
    label = ['1 comp';...
             '2 comp';...
             '3 comp'];
    clabel = cellstr(label);
    lb_pos_x = xc;
    lb_pos_y = [yc(1)-0.04*lc,yc(2)-0.04*lc,yc(3)+0.04*lc];
    for i=1:3
        text(lb_pos_x(i),lb_pos_y(i),clabel(i),...
             'FontSize', 14,...
             'HorizontalAlignment', 'center');
    end

    label_side1 = 'Axisymmetric expansion';
    text((xc(1)+xc(3))/2,(yc(1)+yc(3))/2+0.08*lc,label_side1,...
         'Rotation', -60,...
         'HorizontalAlignment', 'center',...
         'FontSize', 12);
    label_side2 = 'Axisymmetric contraction';
    text((xc(2)+xc(3))/2,(yc(2)+yc(3))/2+0.08*lc,label_side2,...
         'rotation', 60,...
         'HorizontalAlignment', 'center',...
         'FontSize', 12);
    label_side3 = 'Two component';
    text((xc(1)+xc(2))/2,(yc(1)+yc(2))/2-0.04*lc,label_side3,...
         'HorizontalAlignment', 'center',...
         'FontSize', 12);
    label_ps = 'Plane strain';
    lx_ps = abs(x1_ps-x2_ps);
    ly_ps = abs(y1_ps-y2_ps);
    ang = atan2(ly_ps,lx_ps)./pi.*180;
    text((x1_ps+x2_ps)/2,(y1_ps+y2_ps)/2+0.16*lc,label_ps,...
         'Rotation', ang,...
         'HorizontalAlignment', 'center',...
         'FontSize', 12);
end
