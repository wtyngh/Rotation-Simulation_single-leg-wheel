%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is used for observing the difference between leg and wheel
% Especially for the different characteristics on different terrains
%
% Geometry included
% Dynamic condisered
%
% Last advised : 2018/02/01
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opengl info


%% Draw the continuous animation with given conditions 
clear variables; clc;
timer_total = tic;

control_mode = 2;
% 1:torque control, input 'time-torque' trajectory
% 2:Position control, input 'time-theta' trajectory

radius = 0.11;  % leg length (m)
delta_r_initial = 0.045;  % delta leg length (m) [0 , 0.045]
leg_mass = 1 ; % define the mass of the structure (kg)
leg_inertia = 1; % define the inertia of the structure (kg*m^2)

mu_s = 0.9; % define the equivalent static friction constant between the wheel and the ground 
mu_k = 0.8; % define the equivalent dynamic friction constant


mass_force = [0 -(leg_mass*9.8)];
%% Settings 

enable.video = 1;  % switch to 1 to enable video recording
enable.xls_record = 0;   % switch to 1 to write the data to the excel file
enable.time_elapsed_print = 1;  % switch to 1 to show the time elapsed of each iteration
enable.plot_quiver = 1;  % switch to 1 to show the force quiver including mass and reaction force from the ground
enable.plot_required_torque = 1; % switch 1 to show the required_torque
enable.plot_procedure = 1; % switch 1 to plot all the procedure

visualization.force = 0.02; % set the quiver factor for the force vector 
visualization.movement = 25; % set the quiver factor for the movement vector 

%% Inital values 
hip_joint_initial = [0,0.2];  % initail position of the hip joint
theta_initial = 0; % define the intial posture of the leg
theta_end = theta_initial + 2 * pi; % define the fianl posture of the leg

% define how much time the leg is going to run (sec)
t_initial = 0;  % (s)
t_end = 10; 

V_initial = [0 0] ; % (m/s)

% define the resolution of the animation
% More points, higher resolution 
num_of_iterations = 501;


t_array = linspace(t_initial, t_end, num_of_iterations);  % t
% define the gait table
switch control_mode
    case 1  % torque control
        % stay tuned
        torque_array =  2 * ones(size(t_array)); % constant torque output
        
    case 2  % position control
        theta_array = linspace(theta_initial, theta_end, num_of_iterations); % constant omega 
        r_array = 0 * theta_array + delta_r_initial ;  % constant delta_r
        
        omega_array = [diff(theta_array),0];
        alpha_array = [diff(omega_array),0];    
end

t_increment = (t_end - t_initial)/ (num_of_iterations - 1);

%% Define landscape 

x_range = [-0.2, 1.5]; % range of the window
y_range = [-0.2, 0.6];

x_partition_diff = 0.001; % define the resolution of the gound
x_partition = x_range(1):x_partition_diff:x_range(2);  % x_partition

landscape_function_index = 4;  

switch(landscape_function_index)
    case 1   % Rough terrain
        landscape_function = @(x) 0.09*sin(10*x) + x*0.1 ;
        landscape_str = 'rough';
    case 2   % Flat terrain
        landscape_function = @(x) 0 * x   ;
        landscape_str = 'flat';
%     case 3   % Stairs
%         landscape_partition = craete_stair_landscape(x_partition, 6, 8) ;  
%         % (x_partition, stair_level, level_height)
%         landscape_str = 'stairs';
    case 4   % parabolic
        landscape_function = @(x) 0.8 * (x + 0.1).^2  ;  
        landscape_str = 'parabolic';
end

landscape_partition = landscape_function(x_partition);
str_landscape_function = func2str(landscape_function);
landscape_str_full = [landscape_str,' ',str_landscape_function(5:end)];

landscape_partition_diff = [0, diff(landscape_partition)];
% angle between contact point_tangential and horizontal line
landscape_partition_phi = atan2(landscape_partition_diff, x_partition_diff);

landscape_partition_slope = landscape_partition_diff ./ x_partition_diff;


landscape_table = [ x_partition;
                    landscape_partition];

% land_table = [ x_partition;
%                landscape_partition_diff;
%                landscape_partition_slope ];

% clear x_partition landscape_partition landscape_partition_diff;



%% Define dynamic equations 

% E = 10^-8;  % for numeric differential increment

% ===== define trajectory function =====
% landscape function
Fn_phi = @(x) interp1(x_partition, landscape_partition_phi, x,'linear','extrap');
Fn_slpoe = @(x) interp1(x_partition, landscape_partition_phi, x,'linear','extrap');
% theta trajectory function
Fn_theta = @(t) interp1(t_array,theta_array,t,'linear','extrap');
% r trajectory function
Fn_dr = @(t) interp1(t_array,r_array,t,'linear','extrap');

% ====== define geometry =====

% delta theta, initial position at 3/2 pi
Fn_th = @(t) [0,0,-Fn_theta(t)]; 

% from half cercle center to contact point
% reference x is with respect to the contact point
Fn_R = @(x) radius * [sin(Fn_phi(x)),-cos(Fn_phi(x)),0]; 

% from hip to the center half circle_1
Fn_r1 = @(t) Fn_dr(t) * [sin(Fn_theta(t)),cos(Fn_theta(t)),0];
% from hip to the center half circle_2
Fn_r2 = @(t) -Fn_dr(t)*[sin(Fn_theta(t)),cos(Fn_theta(t)),0];  

 
% Position_1 = @(x,t) cross(R(x),th(t)) - ( R(x) + r1(t) );
Fn_Position_1 = @(x,t)  - ( Fn_R(x) + Fn_r1(t) );
Fn_Vel_1 = @(x,t) (Fn_Position_1(x,t+t_increment)-Fn_Position_1(x,t))/t_increment ...
            + cross(Fn_R(x), (Fn_th(t+t_increment)-Fn_th(t))/t_increment )...
            + (Fn_Position_1(x+x_partition_diff,t)-Fn_Position_1(x,t))/x_partition_diff;
       
    
    
Fn_Acc_1 = @(x,t) (Fn_Vel_1(x,t+t_increment)-Fn_Vel_1(x,t))/t_increment ;...
%       + (Fn_Vel_1(x+x_partition_diff,t)-Fn_Vel_1(x,t))/x_partition_diff;

Fn_Position_2 = @(x,t) - ( Fn_R(x) + Fn_r2(t) );
% velocity: differentail of position
% d/dx + d/dt
Fn_Vel_2 = @(x,t) (Fn_Position_2(x,t+t_increment)-Fn_Position_2(x,t))/t_increment ...
          + cross(Fn_R(x), (Fn_th(t+t_increment)-Fn_th(t))/t_increment )...
          + (Fn_Position_2(x+x_partition_diff,t)-Fn_Position_2(x,t))/x_partition_diff;
        
    
    
Fn_Acc_2 = @(x,t) (Fn_Vel_2(x,t+t_increment)-Fn_Vel_2(x,t))/t_increment ;...
%             + (Fn_Vel_2(x+x_partition_diff,t)-Fn_Vel_2(x,t))/x_partition_diff;


%% Video settings
if enable.video == 1
    enable.plot_procedure = 1;
    video_filename = ['T=',num2str(t_end ),'(s)'...
                      ', Theta=',num2str(theta_initial*180/pi),'~',num2str(theta_end*180/pi),'(deg)'...
                      ', dr=',num2str(delta_r_initial),...
                      ', mu_s=',num2str(mu_s),...
                      ', mu_k=',num2str(mu_k),...
                      ', ',landscape_str,'.avi'];
    writerObj = VideoWriter(video_filename);
    writerObj.FrameRate = 1 / t_increment;  % set playing frame rate
    open(writerObj);   
end

%% Plot the landscape and the leg with initial value
if enable.plot_procedure == 1
    figure(1)
    set(gcf,'name','Leg rotaion simulation','Position', [100 100 1500 800]);
else
    enable.plot_quiver = 0;
end
% First trial
% To get the leg_contour for the further contacting calculation
hip_joint = hip_joint_initial;
V_last = V_initial;
leg_contour = def_leg_contour(hip_joint, theta_initial, delta_r_initial);
movement_vector = [0 0];

% initialize data_record
data_record = double.empty(5,0);


%% Main loop start

for loop_iteration = 1:num_of_iterations

    timer_loop = tic;
    
    if enable.plot_procedure == 1
        subplot(5,1,1:4); 
    end
    
    t = t_array(loop_iteration);
    theta = theta_array(loop_iteration);
    omega = omega_array(loop_iteration);
    delta_r = r_array(loop_iteration);
    
    
    % apply the movement
    hip_joint = hip_joint + movement_vector; 
    
    % Return the leg_contour
    leg_contour = def_leg_contour(hip_joint, theta, delta_r);
    
    %% Check overlap and update the hip joint and contact point
    % Geometric constrian check and fix
    
    % Find contact point and the 
    contact_point = find_contact_point(leg_contour , landscape_table , radius);
    
    if ~isempty(contact_point.point_1)
        
        % tangent of normal force point
        land_diff = interp1(x_partition, landscape_partition_diff, contact_point.point_1(1));  
                 
        tangent_direction = [x_partition_diff , land_diff];
        tangent_direction = tangent_direction / norm(tangent_direction);
        revise_direction = [-land_diff , x_partition_diff] ;
        revise_direction = revise_direction / norm(revise_direction);
        
        % visualize the position shifting  by using arrow
        % according to the overlap
%         force_mag = -100 * contact_point.point_1(3);  % scaled parameter for visualization
%         quiver(contact_point.point_1(1),contact_point.point_1(2),...
%                -force_mag * land_diff , force_mag * x_partition_diff,... (-y,x)
%                 'MaxHeadSize',0.5,'color','b');
        
        % plot contacting point
        if enable.plot_procedure == 1
            plot_legend.contact_point_1 = ...
                plot(contact_point.point_1(1),contact_point.point_1(2),'marker','*','MarkerSize',10,'color','b');
            hold on;
        end
        
        
        contact_point_1 = [contact_point.point_1(1) , contact_point.point_1(2)];
        % assume rolling point is contact point 1
        rolling_point.leg_index = 1;
        rolling_point.point = contact_point_1;
        rolling_point.tangent_force_dir = tangent_direction;
        rolling_point.normal_force_dir = revise_direction;
              
%         revise_vector_1 = revise_distance * revise_direction;
        revise_vector_1 = contact_point.point_1(3:4);
        
    else
        contact_point_1 = [];
        revise_vector_1 = [0,0];
    end
    
    
    if ~isempty(contact_point.point_2)
        
        land_diff = interp1(x_partition, landscape_partition_diff, contact_point.point_2(1));       
        
        tangent_direction = [x_partition_diff , land_diff];
        tangent_direction = tangent_direction / norm(tangent_direction);
        revise_direction = [-land_diff , x_partition_diff] ;
        revise_direction = revise_direction / norm(revise_direction);

        % Visualize the position shifting by using arrow
        % according to the overlap        
%         force_mag = -100 * contact_point.point_2(3);  % scaled parameter
%         quiver(contact_point.point_2(1),contact_point.point_2(2),...
%                -force_mag * land_diff, force_mag*x_partition_diff,... (-y,x)
%                 'MaxHeadSize',0.5,'color','r');
        % plot contacting point
        if enable.plot_procedure == 1
            plot_legend.contact_point_2 = ...
                plot(contact_point.point_2(1),contact_point.point_2(2),'marker','*','MarkerSize',10,'color','r');
            hold on;
        end
        
        
        contact_point_2 = [contact_point.point_2(1) , contact_point.point_2(2)];
        revise_vector_2 = contact_point.point_2(3:4);
   
        
        % redecide rolling point w.r.t contact point 2
        if isempty(contact_point_1)
            rolling_point.leg_index = 2;
            rolling_point.point  = contact_point_2;
            rolling_point.tangent_force_dir = tangent_direction;
            rolling_point.normal_force_dir = revise_direction;
            
        elseif contact_point_2(1) > contact_point_1(1)  % Two contact point, the rolling center is the right one
            rolling_point.leg_index = 2;
            rolling_point.point  = contact_point_2;
            rolling_point.tangent_force_dir = tangent_direction;
            rolling_point.normal_force_dir = revise_direction;
            
        end
        
    else
        contact_point_2 = [];
        revise_vector_2 = [0,0];
    end
       
    
    
    if( isempty(contact_point_1) && isempty(contact_point_2) )
        % no contact, should fall
        rolling_point.leg_index = [];
        rolling_point.point = [];
        rolling_point.normal_force_dir = [];
        isContact = false;
    else
        % contact to ground
        isContact = true;
        if enable.plot_procedure == 1
            rolling_point_txt = ['Rolling point = (',num2str(rolling_point.point (1),4),', ',num2str(rolling_point.point (2),4),' )'];
            text(rolling_point.point (1) , rolling_point.point (2) - 0.1, rolling_point_txt,'color', 'k', 'fontsize', 12);

            plot_legend.rolling_point = plot (rolling_point.point (1), rolling_point.point (2),'marker','.','MarkerSize',20,'color','g');
        end
    end
    
    
    % Adjust the hip joint position
    revise_vector = [max(revise_vector_1(1),revise_vector_2(1)) , max(revise_vector_1(2),revise_vector_2(2))];
    hip_joint = hip_joint + revise_vector;
    leg_contour = def_leg_contour(hip_joint, theta, delta_r);
    %% Drawings 
    
    if enable.plot_procedure == 1
        % Draw the landscape and the leg
        plot_legend = plot_landscape_leg(landscape_table,leg_contour);

        title_str = [sprintf('T = %.2f',t), ' (s) , ',...
                    '\Delta \theta = ', sprintf('%.2f',theta*180/pi),' \circ , ',...
                    '\Delta r = ', sprintf('%.1f',delta_r*100),' (cm) , '...
                    '\mu_s = ', sprintf('%.1f',mu_s),...
                    ' , \mu_k = ', sprintf('%.1f',mu_k),...
                    ' , ', landscape_str ];

        title(title_str, 'fontsize',18);
        axis equal;
        axis([x_range y_range]); % acorrding to the given landscape

        V_txt = ['V = (',sprintf('%.2f',V_last(1)),',',...
        sprintf('%.2f',V_last(2)),') , |V| = ',sprintf('%.2f',norm(V_last)),'(m/s)'] ;
        text( x_range(2) - 0.4 , y_range(1) + 0.08 , V_txt ,'color', 'k', 'fontsize', 12);
        
        text( x_range(1) + 0.05 , y_range(2) - 0.05, landscape_str_full,'color', 'k','fontsize', 12)
    end
   
    %% Determin next step : revolution considering slip effect 
    
    % Kinetics considered
    
    % Calculate reaction force
    if(isContact)   % Contact with ground
        switch rolling_point.leg_index

            case 1
                pos = Fn_Position_1(rolling_point.point(1),t);
                position = rolling_point.point + pos(1:2);
                vel = Fn_Vel_1(rolling_point.point(1),t);
                acc = Fn_Acc_1(rolling_point.point(1),t);
            case 2
                pos = Fn_Position_2(rolling_point.point(1),t);
                position = rolling_point.point + pos(1:2);
                vel = Fn_Vel_2(rolling_point.point(1),t);
                acc = Fn_Acc_2(rolling_point.point(1),t);                
        end

        % Assume no slip first, Calculating tangential force
        % if tangential force > F_static_friction (mu_s)  =>  Slipping
            % => tangential force = F_dynamic_friction (mu_k)
        % Else => rolling without slipping
           
        
        % Assume no slip first
        
        % Force equilibrium
        rolling_point.reaction_force = leg_mass * acc(1:2) - mass_force;  % F + W = ma
        
%         
%         if enable.plot_quiver == 1
% 
%             quiver(hip_joint(1),hip_joint(2),...
%             acc(1)*visualization.force , acc(2)*visualization.force,... 
%             'MaxHeadSize',0.5,'color','r', 'LineStyle', '-'); 
%         
%             quiver(hip_joint(1),hip_joint(2),...
%             vel(1)*visualization.force , vel(2)*visualization.force,... 
%             'MaxHeadSize',0.5,'color','b', 'LineStyle', '-'); % brown
%         end
        
 
        rolling_point.tangent_force = dot(rolling_point.reaction_force, rolling_point.tangent_force_dir)*rolling_point.tangent_force_dir;
        rolling_point.tangent_force_dir = rolling_point.tangent_force / norm(rolling_point.tangent_force);
        
        rolling_point.normal_force = dot(rolling_point.reaction_force, rolling_point.normal_force_dir)*rolling_point.normal_force_dir;
        rolling_point.normal_force_dir = rolling_point.normal_force / norm(rolling_point.normal_force);


        % max friction force
        max_static_friction = mu_s * norm(rolling_point.normal_force);
%         max_static_friction_force = max_static_friction * rolling_point.tangent_force_dir ;
        
        % Equivelent R, from contact point to the hip
        rotation_radius_vector = hip_joint - rolling_point.point ; % contact point to the hip

        if( norm(rolling_point.tangent_force) <= max_static_friction )
            % No-slip condition, rolling with respect to the contact point 
            % calculate the total reaction force provided by ground
            % Rolling
            
            isRolling = true;
            if enable.plot_procedure == 1
                text( x_range(1) + 0.05 , y_range(2) - 0.1, 'No slip','color', 'k', 'fontsize', 12);
            end
        else
            % Slip condition
            % calculate the total reaction force provided by ground
            dynamic_friction_force = mu_k * norm(rolling_point.normal_force) * rolling_point.tangent_force_dir;
            rolling_point.tangent_force = dynamic_friction_force;
            rolling_point.reaction_force = rolling_point.normal_force + rolling_point.tangent_force;
            
            
            % transfer the external force to displacement
            isRolling = false;  % not static, considering kinetics
            if enable.plot_procedure == 1
                text( x_range(1) + 0.05 , y_range(2) - 0.1, 'Slipping !','color', 'r','fontsize', 12);
            end 
        end

        
        % visualize the force including mass, reaction normal and reaction tangential
        if enable.plot_quiver == 1

            plot_legend.mass_force = quiver(hip_joint(1),hip_joint(2),...
            mass_force(1) * visualization.force , mass_force(2) * visualization.force,... 
            'MaxHeadSize',0.5,'color','k', 'LineStyle', ':');

            plot_legend.reaction_force = quiver(rolling_point.point(1),rolling_point.point(2),...
            rolling_point.reaction_force(1) * visualization.force , rolling_point.reaction_force(2) * visualization.force,... 
            'MaxHeadSize',0.5,'color','k', 'LineStyle', ':');

            % normal reaction force
            plot_legend.reaction_normal_force = quiver(rolling_point.point(1),rolling_point.point(2),...
            rolling_point.normal_force(1)*visualization.force , rolling_point.normal_force(2)*visualization.force,... 
            'MaxHeadSize',0.5,'color',[0.6350 0.0780 0.1840], 'LineStyle', ':'); % brown

            % tangential reaction force
            plot_legend.reaction_tangent_force = quiver(rolling_point.point(1),rolling_point.point(2),...
            rolling_point.tangent_force(1)*visualization.force , rolling_point.tangent_force(2)*visualization.force,... 
            'MaxHeadSize',0.5,'color',[0.6350 0.0780 0.1840], 'LineStyle', ':'); % brown

        end
  
        % adjust velocity
        if dot(V_last, -rolling_point.normal_force_dir) > 0
            V_last = V_last - dot(V_last, -rolling_point.normal_force_dir)*(-rolling_point.normal_force_dir);
        end

        
    else
        
        % Does not contact to ground, fall.
        rolling_point.reaction_force = [0,0];
        rotation_radius_vector = [0,0];

        isRolling = false;  % not static, considering kinetics
        text(  x_range(1) + 0.05 , y_range(2) - 0.1, 'Falling !','color', 'k', 'fontsize', 12);
    end
    
    total_force = rolling_point.reaction_force + mass_force;
    
    if enable.plot_quiver == 1
        % total force
        plot_legend.total_force = quiver(hip_joint(1),hip_joint(2),...
        total_force(1)*visualization.force , total_force(2)*visualization.force,... 
        'MaxHeadSize',0.5,'color',[0.6350 0.0780 0.1840], 'LineStyle', '--'); % brown
    end
    
    
    acceleration = (total_force) / leg_mass;
    
    required_torque = leg_inertia * alpha_array(loop_iteration)...
                    + cross([rolling_point.reaction_force 0],[rotation_radius_vector 0]);
    required_torque = required_torque(3);
    
    
    if enable.plot_procedure == 1
    
        a_txt = ['a = (',sprintf('%.2f',acceleration(1)),',',...
        sprintf('%.2f',acceleration(2)),') , |a| = ',sprintf('%.2f',norm(acceleration)),'(m/s)'] ;
        text( x_range(2) - 0.4 , y_range(1) + 0.04 , a_txt ,'color', 'k', 'fontsize', 12);
    
    end


    
    % visualize the total force by using arrow
%     if enable.plot_quiver == 1
%         plot_legend.total_force = quiver(hip_joint(1),hip_joint(2),...
%                    visualization.force * rolling_point.reaction_force(1),visualization.force * rolling_point.reaction_force(2),... 
%                     'MaxHeadSize',2,'color','r'); 
%     end

        
    % Determine movement
    
    
    if isRolling == true  % Rolling
        % No-slip condition, rolling with respect to the contact point            
        % rotate clockwise wrt the contact point
        
        new_rotation_radius_vector =  rotation_radius_vector * [cos(-omega) sin(-omega) 
                                                               -sin(-omega) cos(-omega)] ;
        movement_vector = (rolling_point.point + new_rotation_radius_vector) - hip_joint;

        V_now = vel(1:2);
%         V_now = V_last + acceleration * t_increment;
%         movement_vector = position - hip_joint;
    else  
        % not static, considering kinetics
        % additional force convert to acceleration   
        % including falling and slipping
        V_now = V_last + acceleration * t_increment;
        movement_vector = V_last * t_increment + ...
                          0.5 * acceleration * t_increment^2;
    end

    V_last = V_now; 
    
    
    % visualize the hip joint movement by using arrow
    % now hip joint position
    % scaled parameter
    if enable.plot_quiver == 1
        plot_legend.movement = quiver(hip_joint(1),hip_joint(2),...
                       visualization.movement * movement_vector(1),visualization.movement * movement_vector(2),... 
                        'MaxHeadSize',0.5,'color','k');
    end
    
    %% Record data, adjust array size with loop
    data_record(1,loop_iteration) = t;
    data_record(2,loop_iteration) = theta;
    data_record(3,loop_iteration) = hip_joint(1);
    data_record(4,loop_iteration) = hip_joint(2);
    data_record(5,loop_iteration) = required_torque;
    
    if enable.plot_procedure == 1
    
        % plot the trajectory of the hip joint
        plot_legend.hip = plot(data_record(3,:),data_record(4,:),...
                'marker','.','MarkerSize',2,'color',[0.4660   0.6740   0.1880]);

        % plot the legend    
        legend([plot_legend.landscape plot_legend.hip plot_legend.leg_1 plot_legend.leg_2 plot_legend.movement],...
                {'Landscape','Hip joint trajectory','Leg_1','Leg_2','Movement vector'},...
                'FontSize',14);
    end
        
    % write video or refresh drawing
    if enable.video == 1
        videoFrame = getframe(gcf);
        writeVideo(writerObj, videoFrame);
    else
        if enable.plot_procedure == 1
            drawnow;
        end
    end
    hold off;
    
    % print the elapsed time
    if enable.time_elapsed_print == 1
        time_str = [sprintf('%.1f',(loop_iteration/num_of_iterations*100)),'%% , ',...
                    sprintf('Elapsed = %.2f(s)', toc(timer_total))...
                    sprintf(', loop = %.2f(s)\n', toc(timer_loop))];
        fprintf(time_str);
    end
    
    if enable.plot_procedure == 1
        subplot(5,1,5);   
        plot(data_record(1,:),data_record(5,:),'color',[ 0    0.4470    0.7410],'linewidth',1.5);
        hold on;
        plot([0 t_end],[0 0],'--','color',[0.01 0.01 0.01]);
        title(['Minimun torque require = ',sprintf('%.2f',required_torque),' (Nm)']);
        xlabel('time (s)');
        ylabel('Torque (Nm)');
        xlim([t_initial t_end]);
        hold off;
    end
    
    
    
end

%% Calculation
% Calculate min. required work
total_work = trapz( data_record(1,:),data_record(5,:));
total_work_abs = trapz( data_record(1,:) , abs(data_record(5,:)) );

% Calculate the length of hip joint trajectory  
% Calculate integrand from x,y derivatives, and integrate to calculate arc length
hip_joint_trajectory_length =  trapz(hypot(   diff( data_record(3,:) ), diff( data_record(4,:) )  ));   

% Calculate the length of the landscape, where the hip joint has traveled 
traveled_landscape.points = [data_record(3,:) ; 
                            interp1(x_partition, landscape_table(2,:), data_record(3,:) )];
traveled_landscape.length = trapz(hypot( diff(traveled_landscape.points(1,:)) , diff(traveled_landscape.points(2,:)) ));


hip_joint_vs_landscape_length_ratio = hip_joint_trajectory_length / traveled_landscape.length
work_per_landscape_length = total_work_abs / traveled_landscape.length

data_record(6,1) = hip_joint_vs_landscape_length_ratio;
data_record(7,1) = work_per_landscape_length;


%%
if enable.video == 1
    close(writerObj);
    fprintf('video finished\n');
end

if enable.xls_record == 1
    xlsx_tab_str = ['theta=',num2str(theta_end*180/pi),', dr=',num2str(delta_r),', ', landscape_str];
    data_record = data_record'; % switch arrangement from row to column
%     data_col_header = {'T' ,'Theta','Hip joint x','Hip joint y','Min required torque'};

    [xls_status, xls_message] = xlswrite(' 20180117.xlsx',data_record, xlsx_tab_str);
%     [xls_status, xls_message] = writetable(data_table,'table.xlsx','Sheet', xlsx_tab_str);
    if xls_status == 1
        fprintf('xlsx write sucessful\n');
    else
        fprintf('xlsx write error\n');
    end
end

fprintf('Total time = %f sec\n', toc(timer_total));
%%
% figure(2)
% set(gcf,'name','minimun torque require');
% plot(data_record(1,:),data_record(7,:),'linewidth',1.5);
% hold on;
% plot([0 t_end],[0 0],'--','color',[0.01 0.01 0.01]);
% title('Minimun torque require');
% xlabel('time (s)');
% ylabel('Torque (Nm)');