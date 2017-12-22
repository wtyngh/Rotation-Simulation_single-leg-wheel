function assigned_point = get_assigned_position(leg_num, hip_joint, delta_theta, delta_r, assinged_theta)
% return the position [x,y] of the assigned point

    r = 11;  % leg radius
    real_theta = - delta_theta + 3/2*pi ;
    assigned_point = [0,0];
    
    switch (leg_num)
        case 1
            leg_1_center = [hip_joint(1) - delta_r*cos(real_theta), ...
                            hip_joint(2) - delta_r*sin(real_theta)];
                        
            assigned_point = [leg_1_center(1) + r*cos(real_theta + assinged_theta),...
                              leg_1_center(2) + r*sin(real_theta + assinged_theta)];
                        
        case 2
            leg_2_center = [hip_joint(1) + delta_r*cos(real_theta), ...
                            hip_joint(2) + delta_r*sin(real_theta)];
                        
            assigned_point = [leg_2_center(1) + r*cos(real_theta + assinged_theta),...
                              leg_2_center(2) + r*sin(real_theta + assinged_theta)]; 
            
    end
end