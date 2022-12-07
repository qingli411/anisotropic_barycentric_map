function r = setupVectorDirection(varargin)
% setupVectorDirection
%   r = setupVectorDirection() sets up the vector direction map and
%   returns the radius of the direction circle in the map which is used
%   to get the x- and y-coordinates (1 by default).
%   No arguments are required.
%
%   See also vectorDirectionCoord


    nArgs = length(varargin);
    if nArgs > 0
        error('No arguments required.');
    end

    hold on;
    % plot circles
    theta = 0:1:360;
    d2r = pi/180;
    theta = d2r.*theta; % degree to radian
    r = 1;
    circle_color = [0.5, 0.5, 0.5];
    xcircle = r.*cos(theta);
    ycircle = r.*sin(theta);
    plot(xcircle, ycircle,...
         'LineWidth', 2,...
         'LineStyle', '-',...
         'Color', circle_color);
    % plot x-y angles
    nr = 12;
    for i=1:nr
        % draw reference lines
        phi = pi.*2./nr.*i;
        xref = r.*cos(phi);
        yref = r.*sin(phi);
        plot([0 xref], [0, yref],...
             'LineWidth', 1.0,...
             'LineStyle', '--',...
             'Color', circle_color);
        % add label
        if i~=nr && i~=nr/4
            di = i*360/nr;
            if di>180
                di = di-360;
            end
            di_str = sprintf('%i', di);
            dl = 0.15;
            xl = xref+dl*cos(phi);
            yl = yref+dl*sin(phi);
            addLabelDeg(xl, yl, di_str, circle_color, 9);
        end
    end
    % plot z angles
    nz = 6;
    for i=1:nz-1
        % draw referecen lines
        ri = r/nz*i;
        xref = ri.*cos(theta);
        yref = ri.*sin(theta);
        plot(xref, yref,...
             'LineWidth', 1.0,...
             'LineStyle', '--',...
             'Color', circle_color);
        % add label
        di = 90-ri.*90;
        di_str = sprintf('%i',di);
        xl = -ri-0.05;
        yl = -0.05;
        addLabelDeg(xl, yl, di_str, circle_color, 8);
    end

    % plot arrows
    ax = 1.25;
    ay = 1.25;
    plotArrow(0, 0, ax, 0);
    plotArrow(0, 0, 0, ay);

    % add label
    addLabel(ax, -0.1, 'x');
    addLabel(-0.1, ay, 'y');
    addLabel(0.06, 0.08, 'z');

    % turn off axis
    axis off;
    % set aspect ratio to 1:1:1
    daspect([1,1,1]);
end

% local functions
function h = plotArrow(x, y, u, v)
    h = quiver(x, y, u, v,...
               'LineWidth', 1.5,...
               'Color', 'k',...
               'AutoScaleFactor', 1.0,...
               "maxheadsize", 0.2);
end

function tx = addLabelDeg(ax, ay, di_str, color, fontsize)
    tx = text(ax, ay, [ di_str '^\circ'],...
              'Interpreter', 'tex',...
              'FontSize', fontsize,...
              'Color', color,...
              'HorizontalAlignment', 'center');
end

function tx = addLabel(ax, ay, str)
    tx = text(ax, ay, str,...
              'Interpreter', 'tex',...
              'FontSize', 14,...
              'HorizontalAlignment', 'center');
end

