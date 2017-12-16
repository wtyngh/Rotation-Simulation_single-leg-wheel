function boundary = draw_leg(hip_joint, delta_theta, delta_r)

% define initial position of the lag
% couunter-clockwise : + delta_thetatheta
    defined_theta = - delta_theta; 
%     delta_r = 1;

    real_theta = defined_theta + 3/2*pi ;

%     hip_joint = [0,0];
    plot(hip_joint(1),hip_joint(2),'marker','o','MarkerSize',8)

    hold on;
    
    % find the center of the two half circle individually
    leg_1_center = [hip_joint(1) - delta_r*cos(real_theta), ...
                    hip_joint(2) - delta_r*sin(real_theta)];

    leg_2_center = [hip_joint(1) + delta_r*cos(real_theta), ...
                    hip_joint(2) + delta_r*sin(real_theta)];

    boundary_1 = draw_half_circle(leg_1_center(1), leg_1_center(2), real_theta, 1);
    boundary_2 = draw_half_circle(leg_2_center(1), leg_2_center(2), real_theta + pi, 2);
    
    boundary.leg_1 = boundary_1;
    boundary.leg_2 = boundary_2;
end