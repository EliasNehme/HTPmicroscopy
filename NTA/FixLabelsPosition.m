function FixLabelsPosition(axh, azimuth, elevation, shifts_x, shifts_y)
% taken from here: https://www.mathworks.com/matlabcentral/answers/1686239-aligning-axes-labels-in-3d-plots

    % get axis angles
    unitx = [1;0;0];
    unity = [0;1;0];
    unitz = [0;0;1];
    projectedunitx = rotx(elevation) * rotz(-azimuth) * unitx;
    projectedunity = rotx(elevation) * rotz(-azimuth) * unity;
    xlabelangle = atan2d(projectedunitx(3),projectedunitx(1));
    ylabelangle = -(180 - atan2d(projectedunity(3),projectedunity(1)));
    
    % fix labels rotation
    xlabelhandle = axh.XLabel;
    ylabelhandle = axh.YLabel;
    xlabelhandle.Rotation = xlabelangle;
    ylabelhandle.Rotation = ylabelangle;
    
    % get axis limits, calc. mean points and extremes
    xlimits = xlim(axh);
    ylimits = ylim(axh);
    zlimits = zlim(axh);
    xmean = mean(xlimits);
    ymean = mean(ylimits);
    xbottom = xlimits(1);
    % ybottom = ylimits(1);
    ytop = ylimits(end);
    zbottom = zlimits(1);
    
    % set label position
    xlabelhandle.Position = [xmean+shifts_x(1) ytop+shifts_x(2) zbottom+shifts_x(3)];
    ylabelhandle.Position = [xbottom+shifts_y(1) ymean+shifts_y(2) zbottom+shifts_y(3)];
end
