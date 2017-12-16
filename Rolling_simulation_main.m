%% Draw the continuous animation with given conditions
clear variables; clc;

r = 11;
delta_r = 4;

%% Inital values
hip_joint = [0,r];
delta_theta = 0;
leg_index = 1;

delta_theta_end = 3*pi; 

%% Define landscape
x_range = [-25, 150];
y_range = [-5, 50];

x_partition_diff = 0.1;
x_partition = x_range(1):x_partition_diff:x_range(2);  % x_partition


% landscape function
landscape_partition = 5 * sin(0.2 * x_partition) + x_partition*0.2 ;
% landscape_partition = 0 * x_partition;

landscape_partition_diff = diff(landscape_partition);
% landscape_partition = landscape(x_partition);
% landscape_partition = subs(landscape_function,x_partition);

landscape_table(1,:) = x_partition;
landscape_table(2,:) = landscape_partition;
landscape_table(3,1) = 0;
landscape_table(3,2:end) = landscape_partition_diff;

clear x_partition landscape_partition_diff;

 
%% Plot landscape and draw the leg with initial value
figure(1)
set(gcf,'name','leg','Position', [100 100 1500 800]);

% decide the resolution of the animation
delta_theta_increment = pi/100;


% first trial
plot(landscape_table(1,:),landscape_table(2,:),'color','k' ,'linewidth', 2);
hold on;
title_str = ['\Delta \theta = ',num2str(delta_theta*180/pi),' \circ',...
             ' , \Delta r = ',num2str(delta_r),' cm'];
title(title_str);
boundary = draw_leg(hip_joint, delta_theta, delta_r);
hold off;


while delta_theta <= delta_theta_end

% for delta_theta = 0:delta_theta_increment: 3*pi
    %% Check overlap and update the hip joint and contact point
    % Find normal force point
    
    normal_force_point = find_normal_force_point(boundary,landscape_table);
    
    
    if ~isempty(normal_force_point.point_1)
        
        % tangent of normal force point
        land_diff = lookup_table(landscape_table(1,:),landscape_table(3,:),normal_force_point.point_1(1));  
        
        
        % l* sin(theta), where theta is the angle between tangent vector
        % and verticle line
        force_distance = abs(normal_force_point.point_1(3))*( x_partition_diff / sqrt(x_partition_diff^2 + land_diff^2) ); 
        
        
        force_distance = force_distance * 0.5;  % for distance estimation error, more reasonable result
        
        
        
        force_direction = [-land_diff , x_partition_diff];
        force_direction = force_direction/ (sqrt(land_diff^2+x_partition_diff^2));
        
        % visualize the position shifting  by using arrow
        % according to the overlap
        force_mag = -50 * normal_force_point.point_1(3);  %scaled parameter
        quiver(normal_force_point.point_1(1),normal_force_point.point_1(2),...
               -force_mag * land_diff , force_mag * x_partition_diff,... (-y,x)
                'MaxHeadSize',0.5,'color','b');
        hold on;
        plot(normal_force_point.point_1(1),normal_force_point.point_1(2),'marker','*','MarkerSize',10,'color','b');


        % adjust the hip joint
        hip_joint = hip_joint + force_distance * force_direction;
        
        contact_point = [normal_force_point.point_1(1) , normal_force_point.point_1(2)];
        
        contact_point_txt = ['(',num2str(contact_point(1),4),',',num2str(contact_point(2),4),')'];
        text(contact_point(1) , contact_point(2) - 2, contact_point_txt, 'color', 'b');
        
    else
        contact_point = [];
    end
    
    if ~isempty(normal_force_point.point_2)

        land_diff = lookup_table(landscape_table(1,:),landscape_table(3,:),normal_force_point.point_2(1));
        
        % visualize the position shifting by using arrow
        % according to the overlap

        force_distance = abs(normal_force_point.point_2(3))*( x_partition_diff / sqrt(x_partition_diff^2 + land_diff^2) );
        
        
        force_distance = force_distance * 0.5;  % for distance estimation error, more reasonable result
        
        
        force_direction = [-land_diff , x_partition_diff];
        force_direction = force_direction / (sqrt(land_diff^2+x_partition_diff^2));
        
        
        force_mag = -50 * normal_force_point.point_2(3);  % scaled parameter
        quiver(normal_force_point.point_2(1),normal_force_point.point_2(2),...
               -force_mag * land_diff, force_mag*x_partition_diff,... (-y,x)
                'MaxHeadSize',0.5,'color','r');            
        hold on;
        plot(normal_force_point.point_2(1),normal_force_point.point_2(2),'marker','*','MarkerSize',10,'color','r');

        
        % adjust the hip joint
        hip_joint = hip_joint + force_distance * force_direction;
        
        if isempty(contact_point) 
            contact_point = [normal_force_point.point_2(1) , normal_force_point.point_2(2)];
            
            contact_point_txt = ['(',num2str(contact_point(1),4),',',num2str(contact_point(2),4),')'];
            text(contact_point(1) , contact_point(2) - 2, contact_point_txt,'color', 'r');
            
        elseif normal_force_point.point_2(1) > contact_point(1)  % the point is on the righter side
            contact_point = [normal_force_point.point_2(1) , normal_force_point.point_2(2)];
            
            contact_point_txt = ['(',num2str(contact_point(1),4),',',num2str(contact_point(2),4),')'];
            text(contact_point(1) , contact_point(2) - 2, contact_point_txt,'color', 'r');
        end
        
    end
    
    
    if ~isempty(contact_point)
        plot (contact_point(1),contact_point(2),'marker','.','MarkerSize',20,'color','c')
        hold on;
    end
    
    
    %% Drawings 
    % draw the landscape
    plot(landscape_table(1,:),landscape_table(2,:),'color','k' ,'linewidth', 2);
    hold on;
    
    
    title_str = ['\Delta \theta = ',num2str(delta_theta*180/pi),' \circ',...
                 ' , \Delta r = ',num2str(delta_r),' cm'];
    title(title_str);
    
    % draw the leg
    % return the contour of the leg boundary
    boundary = draw_leg(hip_joint, delta_theta, delta_r);
    axis equal;
    axis([x_range y_range]);

%     plot(assigned_point(1),assigned_point(2),'marker','*','MarkerSize',10)
%     drawnow;
    

    
    %% Determin next step : revolution using no-slip assumption 
    
%     hip_joint(1) = hip_joint(1) + delta_theta_increment*r ;
    if(~isempty(contact_point))
        % Contact with ground, rolling with respect to the contact point
        
        
        rotation_radius_vector = hip_joint - contact_point;
        rotation_radius_vector_length = sqrt(rotation_radius_vector(1)^2 + rotation_radius_vector(2)^2);
        % Revolute_direction = (y,-x) of rotation_radius_vector(x,y)
        
        new_rotation_radius_vector =  rotation_radius_vector * [cos(-delta_theta_increment) sin(-delta_theta_increment) 
                                                               -sin(-delta_theta_increment) cos(-delta_theta_increment)] ;
        movement_vector = (contact_point + new_rotation_radius_vector) - hip_joint;

%         revolute_direction = [rotation_radius_vector(2),-rotation_radius_vector(1)];
%         revolute_direction = revolute_direction * [cos(-delta_theta_increment/2) sin(-delta_theta_increment/2) 
%                                                    -sin(-delta_theta_increment/2) cos(-delta_theta_increment/2)] ;
%         
%         revolute_direction = revolute_direction / sqrt(revolute_direction(1)^2 + revolute_direction(1)^2);
%         movement_vector = ( delta_theta_increment * rotation_radius_vector_length ) * revolute_direction;
        
        
        delta_theta = delta_theta + delta_theta_increment;
        
    else
        % Does not contact to ground, fall and does not roll.
        
        movement_vector = [0 -0.05];
    end
    
    % visualize the hip joint movement by using arrow
    % now hip joint position
    % scaled parameter
    quiver(hip_joint(1),hip_joint(2),...
           50*movement_vector(1),50*movement_vector(2),... 
            'MaxHeadSize',0.5,'color','k'); 

    hip_joint = hip_joint + movement_vector;
    



%     assigned_point = get_assigned_position(leg_index, hip_joint, delta_theta, delta_r, delta_theta);

    
    drawnow;
    hold off;

end



