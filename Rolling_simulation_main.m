%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is used for observing the difference between leg and wheel
% Especially for the different characteristics on different terrains
%
% Geometry included
%
% Last advised : 2018/01/11
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opengl info


%% Draw the continuous animation with given conditions
clear variables; clc;

r = 0.11;  % leg length (m)
delta_r = 0;  % delta leg length (m)
mass = 0.1 ; % define the mass of the structure (kg)

static_friction_const = 0.9; % define the equivalent static friction constant between the wheel and the ground 
mass_vector = [0 -(mass*9.8)];

enable_video = 0;  % switch to 1 to enable video recording
enable_xls_record = 0;   % switch to 1 to write the data to the excel file
%% Inital values
hip_joint = [0,0.15];
delta_theta = 0;
leg_index = 1;

% define how much the leg is going to run
delta_theta_end = 2 * pi; 

%% Define landscape
x_range = [-0.2, 1.5];
y_range = [-0.15, 0.6];

x_partition_diff = 0.01;
x_partition = x_range(1):x_partition_diff:x_range(2);  % x_partition

landscape_function = 1;

switch(landscape_function)
% landscape function
    case 1   % Rough terrain
        landscape_partition = 0.05 * sin(10* x_partition) + x_partition*0.8 ;
        landscape_str = 'rough';
    case 2   % Flat terrain
        landscape_partition = 0 * x_partition  ;
        landscape_str = 'flat';
    case 3   % Stairs
        landscape_partition = craete_stair_landscape(x_partition, 6, 8) ;  
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

%% Plot the landscape and the leg with initial value
figure(1)
set(gcf,'name','Leg rotaion simulation','Position', [100 100 1500 800]);

% decide the resolution of the animation; 
% smaller increment, higher resolution
delta_theta_increment = pi/200;

% First trial
% To get the leg_contour for the further contacting calculation
leg_contour = def_leg_contour(hip_joint, delta_theta, delta_r);
plot_legend = plot_landscape_leg(landscape_table,leg_contour);
% plot(landscape_table(1,:),landscape_table(2,:),'color','k' ,'linewidth', 2);

title_str = ['\Delta \theta = ',num2str(delta_theta*180/pi),' \circ',...
             ' , \Delta r = ',num2str(delta_r),' cm , ', landscape_str];
title(title_str);
axis equal;
axis([x_range y_range]); % acorrding to the given landscape

hold off;

%% Main loop start

loop_iteration = 1;
while delta_theta <= delta_theta_end
    %% Check overlap and update the hip joint and contact point
    
    % Find normal force point
    normal_force_point = find_normal_force_point(leg_contour,landscape_table);
    
    if ~isempty(normal_force_point.point_1)
        
        % tangent of normal force point
        land_diff = lookup_table(landscape_table(1,:),landscape_table(3,:),normal_force_point.point_1(1));  
                
        % l* sin(theta), where theta is the angle between tangent vector and verticle line
        force_distance = abs(normal_force_point.point_1(3))*( x_partition_diff / sqrt(x_partition_diff^2 + land_diff^2) ); 
        
        % The steeper slope, the smaller value, range(0,1]
        force_distance = force_distance * (0.5);  % for distance estimation error, more reasonable result
        
        force_direction = [-land_diff , x_partition_diff];
        % normalize
        force_direction = force_direction/ (sqrt(land_diff^2+x_partition_diff^2));
        
        % visualize the position shifting  by using arrow
        % according to the overlap
        force_mag = -100 * normal_force_point.point_1(3);  % scaled parameter for visualization
        quiver(normal_force_point.point_1(1),normal_force_point.point_1(2),...
               -force_mag * land_diff , force_mag * x_partition_diff,... (-y,x)
                'MaxHeadSize',0.5,'color','b');
        hold on;
        
        % plot contacting point
        plot_legend.contact_point_1 = ...
            plot(normal_force_point.point_1(1),normal_force_point.point_1(2),'marker','*','MarkerSize',10,'color','b');
        
        contact_point_1 = [normal_force_point.point_1(1) , normal_force_point.point_1(2)];
        rolling_point.point = contact_point_1;
        rolling_point.normal_force_dir = force_direction;
        
        % adjust the hip joint
        hip_joint = hip_joint + force_distance * force_direction;
        
    else
        contact_point_1 = [];
        rolling_point.point = [];
        rolling_point.normal_force_dir = [];
    end
    
    
    if ~isempty(normal_force_point.point_2)

        land_diff = lookup_table(landscape_table(1,:),landscape_table(3,:),normal_force_point.point_2(1));       
        % Visualize the position shifting by using arrow
        % according to the overlap

        force_distance = abs(normal_force_point.point_2(3))*( x_partition_diff / norm([x_partition_diff, land_diff]) );
        
        
%         normal_force_point_slope_exp = exp(-abs(land_diff / x_partition_diff))
        
        force_distance = force_distance *(0.5);  % for distance estimation error, more reasonable result
        
        
        force_direction = [-land_diff , x_partition_diff];
        force_direction = force_direction / (norm([land_diff, x_partition_diff]));
        
        
        force_mag = -100 * normal_force_point.point_2(3);  % scaled parameter
        quiver(normal_force_point.point_2(1),normal_force_point.point_2(2),...
               -force_mag * land_diff, force_mag*x_partition_diff,... (-y,x)
                'MaxHeadSize',0.5,'color','r');            
        hold on;
        % plot contacting point
        plot_legend.contact_point_2 = ...
            plot(normal_force_point.point_2(1),normal_force_point.point_2(2),'marker','*','MarkerSize',10,'color','r');

        
        contact_point_2 = [normal_force_point.point_2(1) , normal_force_point.point_2(2)];
        rolling_point.point = contact_point_2;
        
        % Adjust the hip joint
        hip_joint = hip_joint + force_distance * force_direction;
        
        % redecide rolling center
        if isempty(contact_point_1) 
            rolling_point.point  = contact_point_2;
            rolling_point.normal_force_dir = force_direction;
            
        elseif contact_point_2(1) > contact_point_1(1)  % Two contact point, the rolling center is the right one
            rolling_point.point  = contact_point_2;
            rolling_point.normal_force_dir = force_direction;
            
        end
        
    else
        contact_point_2 = [];
    end
       
    
    if( isempty(contact_point_1) && isempty(contact_point_2) )
        rolling_point.point  = [];
        rolling_point.normal_force_dir = [];
    else
        rolling_point_txt = ['Rolling point = (',num2str(rolling_point.point (1),4),', ',num2str(rolling_point.point (2),4),' )'];
        text(rolling_point.point (1) , rolling_point.point (2) - 0.05, rolling_point_txt,'color', 'k');

        plot_legend.rolling_point = plot (rolling_point.point (1), rolling_point.point (2),'marker','.','MarkerSize',20,'color','g');
    end
    
    
    
    % Adjust array size with loop
    hip_joint_record(1,loop_iteration) = delta_theta;
    hip_joint_record(2,loop_iteration) = hip_joint(1);
    hip_joint_record(3,loop_iteration) = hip_joint(2);
    
    
    %% Drawings 

    % Return the leg_contour
    leg_contour = def_leg_contour(hip_joint, delta_theta, delta_r);
    
    % Draw the landscape and the leg
    plot_legend = plot_landscape_leg(landscape_table,leg_contour);
    % plot(landscape_table(1,:),landscape_table(2,:),'color','k' ,'linewidth', 2);

    title_str = ['\Delta \theta = ',num2str(delta_theta*180/pi),' \circ',...
             ' , \Delta r = ',num2str(delta_r),' cm , ', landscape_str];
    title(title_str);
    axis equal;
    axis([x_range y_range]); % acorrding to the given landscape

    
    % plot the trajectory of the hip joint
    plot_legend.hip = plot(hip_joint_record(2,:),hip_joint_record(3,:),...
                'marker','.','MarkerSize',2,'color',[0.4660  0.6740  0.1880]);
    


    %% Determin next step : revolution considering slip effect 
    
    if(~isempty(rolling_point.point))   % Contact with ground
        
        % considering friction effect
        
        rolling_point.normal_force = ...
            dot(-mass_vector, rolling_point.normal_force_dir) * rolling_point.normal_force_dir;
        
        rolling_point.tangent_force = (-mass_vector) - rolling_point.normal_force;
        
        % max friction force
        max_static_friction = static_friction_const * norm(rolling_point.normal_force);
        rolling_point.tangent_force_dir = (rolling_point.tangent_force) / norm(rolling_point.tangent_force);
        max_static_friction_force = max_static_friction * rolling_point.tangent_force_dir ;
        
        
%         quiver(rolling_point.point(1),rolling_point.point(2),...
%         -mass_vector(1),-mass_vector(2),... 
%         'MaxHeadSize',0.5,'color','k');
%         
%         quiver(rolling_point.point(1),rolling_point.point(2),...
%         rolling_point.normal_force(1),rolling_point.normal_force(2),... 
%         'MaxHeadSize',0.5,'color','m');
%         
%         quiver(rolling_point.point(1),rolling_point.point(2),...
%         rolling_point.tangent_force(1),rolling_point.tangent_force(2),... 
%         'MaxHeadSize',0.5,'color','g');
%     
%         quiver(rolling_point.point(1),rolling_point.point(2),...
%         max_static_friction_force(1), max_static_friction_force(2),... 
%         'MaxHeadSize',1,'color','y');
        
        
        
        
        
        
        if( norm(rolling_point.tangent_force) <= max_static_friction )
            % No-slip condition, rolling with respect to the contact point
              
            rotation_radius_vector = hip_joint - rolling_point.point ; % point to the hip
            rotation_radius_vector_length = norm(rotation_radius_vector);

            new_rotation_radius_vector =  rotation_radius_vector * [cos(-delta_theta_increment) sin(-delta_theta_increment) 
                                                                   -sin(-delta_theta_increment) cos(-delta_theta_increment)] ;
            movement_vector = (rolling_point.point  + new_rotation_radius_vector) - hip_joint;


            delta_theta = delta_theta + delta_theta_increment; %increment delta theta
        else
            % Slip condition, additional force convert to acceleration
            
            
            movement_vector = ...
            0.02*(max_static_friction - norm(rolling_point.tangent_force)) * rolling_point.tangent_force_dir;
            
            delta_theta = delta_theta + delta_theta_increment; %increment delta theta
            
            text(rolling_point.point (1) , rolling_point.point (2) + 0.2, 'SLIP !','color', 'r');
            
        end
        
    else
        % Does not contact to ground, fall and does not roll.
        
%         movement_vector = [0 -0.05];
        movement_vector = [0 -(delta_theta_increment * r)]; 
        %synchronize the falling speed with respect to the forwarding speed
    end
    
    
    
    
    
    
    % visualize the hip joint movement by using arrow
    % now hip joint position
    % scaled parameter
    plot_legend.vel_vector = quiver(hip_joint(1),hip_joint(2),...
                   20 * movement_vector(1),20 * movement_vector(2),... 
                    'MaxHeadSize',0.5,'color','k'); 

    hip_joint = hip_joint + movement_vector;
    
    loop_iteration = loop_iteration + 1;

        
    legend([plot_legend.landscape plot_legend.hip plot_legend.leg_1 plot_legend.leg_2 plot_legend.vel_vector],...
            'Landscape','Hip joint trajectory','Leg_1','Leg_2','Velocity vector');

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

