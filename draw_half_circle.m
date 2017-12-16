function boundary = draw_half_circle(center_x, center_y, initial_angle, half_leg_num)
%% plot a circle


% th = 0:pi/50:2*pi;
%x = 0;
%y = 0;
%initial_angle = 0;

%%
    r = 11;  % leg radius
    th = linspace(initial_angle, initial_angle + pi, 100);
    xunit = r * cos(th) + center_x;
    yunit = r * sin(th) + center_y;

    switch(half_leg_num)
        case 1 
            plot(xunit, yunit, 'color','b' ,'linewidth', 1.5);
            hold on;
            plot([center_x,xunit(1)], [center_y,yunit(1)], 'color','b','linewidth', 1.5);
            plot(center_x,center_y, 'marker','x','MarkerEdgeColor','b', 'MarkerSize',6);
        case 2
            plot(xunit, yunit, 'color','r' ,'linewidth', 1.5);
            hold on;
            plot([center_x,xunit(1)], [center_y,yunit(1)], 'color','r','linewidth', 1.5);
            plot(center_x,center_y, 'marker','x','MarkerEdgeColor','r', 'MarkerSize',6);

    end
    
    boundary(1,:) = xunit;
    boundary(2,:) = yunit;
end