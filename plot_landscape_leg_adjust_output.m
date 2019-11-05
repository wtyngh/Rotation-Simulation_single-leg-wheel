function plot_legend = plot_landscape_leg_adjust_output(landscape_table,leg_contour)
    %% plot landscape
    plot_legend.landscape = plot(landscape_table(1,:),landscape_table(2,:),'color','k' ,'linewidth', 2 );
    hold on;
    
    %% plot hip joint
    plot(leg_contour.center(1),leg_contour.center(2),'o','MarkerEdgeColor','k',...
        'MarkerFaceColor',[.49 1 .63], 'MarkerSize',5);
    
    %% plot leg_1
    plot_legend.leg_1 = ...
    plot(leg_contour.leg_1.contour.x, leg_contour.leg_1.contour.y, 'color','b' ,'linewidth', 1.2);
    
    % plot from center to the first point of the contour
    plot([leg_contour.leg_1.center(1), leg_contour.leg_1.contour.x(1)],...
         [leg_contour.leg_1.center(2), leg_contour.leg_1.contour.y(1)], 'color','b','linewidth', 1.2);
     
    % plot the center of the half circle
    plot(leg_contour.leg_1.center(1), leg_contour.leg_1.center(2),... 
        'marker','x','color','b', 'MarkerSize',6);
    
    %% plot leg_2
    plot_legend.leg_2 = ...
    plot(leg_contour.leg_2.contour.x, leg_contour.leg_2.contour.y, 'color','r' ,'linewidth', 1.2);
    
    % plot from center to the first point of the contour
    plot([leg_contour.leg_2.center(1), leg_contour.leg_2.contour.x(1)],...
         [leg_contour.leg_2.center(2), leg_contour.leg_2.contour.y(1)], 'color','r','linewidth', 1.2);
    
    % plot the center of the half circle
    plot(leg_contour.leg_2.center(1), leg_contour.leg_2.center(2),...   
        'marker','x','color','r', 'MarkerSize',6);

    
    %% show the position of the hip joint
    
    hip_joint_txt = ['Hip joint = (',num2str(leg_contour.center(1),4),', ',num2str(leg_contour.center(2),4),' )'];
%     text(leg_contour.center(1) , leg_contour.center(2) + 0.2, hip_joint_txt,...
%         'color',[0 0.5 0], 'fontsize', 12);
    
end