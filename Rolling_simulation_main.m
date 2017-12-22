%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is used for observing the difference between leg and wheel
% Especially for the different characteristics on different terrains
%
% Geometry included
%
% Last advised : 2017/12/23
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Draw the continuous animation with given conditions
clear variables; clc;

r = 11;  % leg length
delta_r = 4.5;  % delta leg length

enable_video = 0;  % switch to 1 to enable video recording
enable_xls_record = 0;   % switch to 1 to write the data to the excel file
%% Inital values
hip_joint = [0,r];
delta_theta = 0;
leg_index = 1;

% define how much the leg is going to run
delta_theta_end = 2 * pi; 

%% Define landscape
x_range = [-25, 150];
y_range = [-20, 70];

x_partition_diff = 0.1;
x_partition = x_range(1):x_partition_diff:x_range(2);  % x_partition

landscape_function = 3;

switch(landscape_function)
% landscape function
    case 1   % Rough terrain
        landscape_partition = 5 * sin(0.2 * x_partition) + x_partition*0.2 - 5;
        landscape_str = 'rough';
    case 2   % Flat terrain
        landscape_partition = 0 * x_partition - 5 ;
        landscape_str = 'flat';
    case 3   % Stairs
        landscape_partition = craete_stair_landscape(x_partition, 6, 8) - 10;  
        % (x_partition, stair_level, level_height)
        landscape_str = 'stairs';
end

landscape_partition_diff = diff(landscape_partition);

landscape_table(1,:) = x_partition;
landscape_table(2,:) = landscape_partition;

% first dirivative value
landscape_table(3,1) = 0;
landscape_table(3,2:end) = landscape_partition_diff;

clear x_partition landscape_partition_diff;

%% Video settings
if enable_video == 1
    video_filename = ['Theta=',num2str(delta_theta_end*180/pi),...
                      ', dr=',num2str(delta_r),', ', landscape_str,'.avi'];
    writerObj = VideoWriter(video_filename);
    writerObj.FrameRate = 30;  % set playing frame rate
    open(writerObj);   
end

%% Plot landscape and draw the leg with initial value
figure(1)
set(gcf,'name','Leg rotaion simulation','Position', [100 100 1500 800]);

% decide the resolution of the animation; 
% smaller increment, higher resolution
delta_theta_increment = pi/200;

% First trial
% To get the boundary for the further contacting calculation
plot(landscape_table(1,:),landscape_table(2,:),'color','k' ,'linewidth', 2);
hold on;
title_str = ['\Delta \theta = ',num2str(delta_theta*180/pi),' \circ',...
             ' , \Delta r = ',num2str(delta_r),' cm , ', landscape_str];
title(title_str);
boundary = draw_leg(hip_joint, delta_theta, delta_r);
hold off;

%% Main loop start

loop_iteration = 1;
while delta_theta <= delta_theta_end
    %% Check overlap and update the hip joint and contact point
    
    % Find normal force point
    normal_force_point = find_normal_force_point(boundary,landscape_table);
    
    if ~isempty(normal_force_point.point_1)
        
        % tangent of normal force point
        land_diff = lookup_table(landscape_table(1,:),landscape_table(3,:),normal_force_point.point_1(1));  
                
        % l* sin(theta), where theta is the angle between tangent vector
        % and verticle line
        force_distance = abs(normal_force_point.point_1(3))*( x_partition_diff / sqrt(x_partition_diff^2 + land_diff^2) ); 
        
        % The steeper slope, the smaller value, range(0,1]
%         normal_force_point_slope_exp = exp(-abs(land_diff / x_partition_diff));
        force_distance = force_distance * (0.5);  % for distance estimation error, more reasonable result
        
        force_direction = [-land_diff , x_partition_diff];
        % normalize
        force_direction = force_direction/ (sqrt(land_diff^2+x_partition_diff^2));
        
        % visualize the position shifting  by using arrow
        % according to the overlap
        force_mag = -50 * normal_force_point.point_1(3);  % scaled parameter for visualization
        quiver(normal_force_point.point_1(1),normal_force_point.point_1(2),...
               -force_mag * land_diff , force_mag * x_partition_diff,... (-y,x)
                'MaxHeadSize',0.5,'color','b');
        hold on;
        plot(normal_force_point.point_1(1),normal_force_point.point_1(2),'marker','*','MarkerSize',10,'color','b');

        % adjust the hip joint
        hip_joint = hip_joint + force_distance * force_direction;
        
        contact_point = [normal_force_point.point_1(1) , normal_force_point.point_1(2)];
        
        contact_point_txt = ['contact point = (',num2str(contact_point(1),4),', ',num2str(contact_point(2),4),' )'];
        text(contact_point(1) , contact_point(2) - 10, contact_point_txt, 'color', 'b');
        
    else
        contact_point = [];
    end
    
    if ~isempty(normal_force_point.point_2)

        land_diff = lookup_table(landscape_table(1,:),landscape_table(3,:),normal_force_point.point_2(1));
        
        % Visualize the position shifting by using arrow
        % according to the overlap

        force_distance = abs(normal_force_point.point_2(3))*( x_partition_diff / sqrt(x_partition_diff^2 + land_diff^2) );
        
        
%         normal_force_point_slope_exp = exp(-abs(land_diff / x_partition_diff))
        
        force_distance = force_distance *(0.5);  % for distance estimation error, more reasonable result
        
        
        force_direction = [-land_diff , x_partition_diff];
        force_direction = force_direction / (sqrt(land_diff^2+x_partition_diff^2));
        
        
        force_mag = -50 * normal_force_point.point_2(3);  % scaled parameter
        quiver(normal_force_point.point_2(1),normal_force_point.point_2(2),...
               -force_mag * land_diff, force_mag*x_partition_diff,... (-y,x)
                'MaxHeadSize',0.5,'color','r');            
        hold on;
        plot(normal_force_point.point_2(1),normal_force_point.point_2(2),'marker','*','MarkerSize',10,'color','r');

        
        % Adjust the hip joint
        hip_joint = hip_joint + force_distance * force_direction;
        
        if isempty(contact_point) 
            contact_point = [normal_force_point.point_2(1) , normal_force_point.point_2(2)];
            
            contact_point_txt = ['contact point = (',num2str(contact_point(1),4),', ',num2str(contact_point(2),4),' )'];
            text(contact_point(1) , contact_point(2) - 10, contact_point_txt,'color', 'r');
            
        elseif normal_force_point.point_2(1) > contact_point(1)  % the point is on the righter side
            contact_point = [normal_force_point.point_2(1) , normal_force_point.point_2(2)];
            
            contact_point_txt = ['contact point = (',num2str(contact_point(1),4),', ',num2str(contact_point(2),4),' )'];
            text(contact_point(1) , contact_point(2) - 10, contact_point_txt,'color', 'r');
        end
        
    end
    
    
    if ~isempty(contact_point)
        p_contact_point = plot (contact_point(1),contact_point(2),'marker','.','MarkerSize',20,'color','c');
        hold on;
    end
    
    % Adjust array size with loop
    hip_joint_record(1,loop_iteration) = delta_theta;
    hip_joint_record(2,loop_iteration) = hip_joint(1);
    hip_joint_record(3,loop_iteration) = hip_joint(2);
    
    
    %% Drawings 
    % Draw the landscape
    p_landscape = plot(landscape_table(1,:),landscape_table(2,:),'color','k' ,'linewidth', 2 );
    hold on;
    
    p_hip_joint = plot(hip_joint_record(2,:),hip_joint_record(3,:),...
                'marker','.','MarkerSize',2,'color',[0.4660  0.6740  0.1880]);
    
    
    title_str = ['\Delta \theta = ',num2str(delta_theta*180/pi),' \circ',...
                 ' , \Delta r = ',num2str(delta_r),' cm , ', landscape_str];
    title(title_str);
    
    % Draw the leg
    % Return the contour of the leg boundary
    boundary = draw_leg(hip_joint, delta_theta, delta_r);
    axis equal;
    axis([x_range y_range]);
    
    hip_joint_txt = ['Hip joint = (',num2str(hip_joint(1),4),', ',num2str(hip_joint(2),4),' )'];
    text(hip_joint(1) , hip_joint(2) + 15, hip_joint_txt,'color', 'k');
    
%     plot(assigned_point(1),assigned_point(2),'marker','*','MarkerSize',10)


    %% Determin next step : revolution using no-slip assumption 
    
    if(~isempty(contact_point))
        % Contact with ground, rolling with respect to the contact point
              
        rotation_radius_vector = hip_joint - contact_point;
        rotation_radius_vector_length = sqrt(rotation_radius_vector(1)^2 + rotation_radius_vector(2)^2);
        
        new_rotation_radius_vector =  rotation_radius_vector * [cos(-delta_theta_increment) sin(-delta_theta_increment) 
                                                               -sin(-delta_theta_increment) cos(-delta_theta_increment)] ;
        movement_vector = (contact_point + new_rotation_radius_vector) - hip_joint;
        
        
        delta_theta = delta_theta + delta_theta_increment; %increment delta theta
        
    else
        % Does not contact to ground, fall and does not roll.
        
%         movement_vector = [0 -0.05];
        movement_vector = [0 -(delta_theta_increment * r)];
    end
    
    % visualize the hip joint movement by using arrow
    % now hip joint position
    % scaled parameter
    p_velocity_vector = quiver(hip_joint(1),hip_joint(2),...
                   20*movement_vector(1),20*movement_vector(2),... 
                    'MaxHeadSize',0.5,'color','k'); 

    hip_joint = hip_joint + movement_vector;
    
    loop_iteration = loop_iteration + 1;
    
%     assigned_point = get_assigned_position(leg_index, hip_joint, delta_theta, delta_r, delta_theta);

        
    legend([p_landscape p_hip_joint p_velocity_vector ],...
            'Landscape','Hip joint trajectory','Velocity vector');

    drawnow;
    hold off;
    

    
    if enable_video == 1
        videoFrame = getframe(gcf);
        writeVideo(writerObj, videoFrame);
    end

end
if enable_video == 1
    close(writerObj);
    fprintf('video finished\n');
end

if enable_xls_record == 1
    xlsx_tab_str = ['theta=',num2str(delta_theta_end*180/pi),', dr=',num2str(delta_r),', ', landscape_str];
    hip_joint_record = hip_joint_record';

    xlswrite('Rolling comparison X-dir.xlsx',hip_joint_record, xlsx_tab_str);
    fprintf('xlsx finished\n');
end

