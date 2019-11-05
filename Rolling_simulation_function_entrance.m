%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is used for observing the difference between leg and wheel
% Especially for the different characteristics on different terrains
% For exploring data points
%
% Geometry included
% Dynamic condisered
%
% Last advised : 2018/05/17
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opengl info


%% Draw the continuous animation with given conditions 
clear variables; clc;
timer_total = tic;

trajectory_mode = 2;  
% 1:constant omega, constant dr = 0, wheel mode 
% 2:constant omega, constant dr = 0.045, legged mode 
% 3:assigned trajectory, CPG trajectory is used here, Trot
% 4:assigned trajectory, CPG trajectory is used here, Walk

leg_mass = 5 ; % define the mass of the structure (kg)
% define the inertia of the structure (kg*m^2) predict by solidworks
% wheel:0.1118(kg*m^2)
% leg:0.02211(kg*m^2)
% motor inertia is small enough to be neglected

mu_s = 0.9; % define the equivalent static friction constant between the wheel and the ground 
mu_k = 0.8; % define the equivalent dynamic friction constant

amp_var_size = 1;
record_data_row_index = 1;

%% define how much time the leg is going to run (sec)
t_initial = 0;  % (s)
t_end = 5; 
% define the resolution of the animation
% More points, higher resolution 
num_of_iterations = 1001;
t_array = linspace(t_initial, t_end, num_of_iterations);  % t

% initialize data_record
record_data = double.empty(0,12);


%% Define landscape 

landscpae_var.x_range = [-0.5, 4.5]; % range of the window
landscpae_var.y_range = [-0.25, 0.8];
landscpae_var.x_partition_diff = 0.001; % define the resolution of the gound

% amp * sin(freq * x) + bias

% for amp_var_ind = 1:1:amp_var_size
% for amp_var_ind = 1:1:amp_var_size



landscpae_var.amp = 0.1;
landscpae_var.freq = 0.2;
landscpae_var.bias = 0;

landscpae_var.mu_s = mu_s;
landscpae_var.mu_k = mu_k;
    

%% Inital values 

hip_joint_initial = [0,0.2];  % initail position of the hip joint    
theta_initial_assigned_deg_array = 0;
% theta_initial_assigned_deg_array = 0:10:179;  % 18 points
size_different_initial = size(theta_initial_assigned_deg_array,2);

for theta_initial_assigned_deg = theta_initial_assigned_deg_array(1:end)

    theta_initial_assigned_rad = theta_initial_assigned_deg/180*pi;

    % define the gait table                    
    switch trajectory_mode

        case {1,2} % constant omega, constant dr
            forward_vel_set = 0.4; % define forward velocity
            if trajectory_mode == 1 % wheel mode
                delta_r_initial = 0;
                forward_dis = forward_vel_set / 0.11*(t_end-t_initial);
                theta_end = theta_initial_assigned_rad + forward_dis; %(V/r)*t=w*t
                trajectory.leg_inertia = 0.1118; % wheel
            else  % trajectory_mode == 2 % legged mode
                delta_r_initial = 0.045;
                forward_dis = forward_vel_set / 0.155*(t_end-t_initial);
                theta_end = theta_initial_assigned_rad + forward_dis; %(V/r)*t=w*t
                trajectory.leg_inertia = 0.02211; %leg
            end
            % define the gait table
            theta_array_full_shifted = linspace(theta_initial_assigned_rad, theta_end, num_of_iterations); % constant omega 
    %                 theta_array = rem(theta_array_full, 2*pi); % project to [0,2*pi)
            r_array_shifted = 0 * theta_array_full_shifted + delta_r_initial ;  % constant delta_r
            input_trajectory_data_filename = ['const w=', num2str(forward_dis/(t_end-t_initial),'%.1f'),'[rad/s]'];
            trajectory.name = input_trajectory_data_filename;

        case {3,4}  % input assigned trajectory, CPG trajectory is used here
            % Load trajectory data                    
            input_trajectory_data_filename = 'CPG trajectory';
            if trajectory_mode == 3
                input_xlsx_tab_str = 'Trot, V=400';
            else  % trajectory_mode == 4
                input_xlsx_tab_str = 'Walk, V=400';     
            end
            trajectory.name = [input_trajectory_data_filename,', ',input_xlsx_tab_str,'[mm/s]'];
            input_trajectory_data = xlsread([input_trajectory_data_filename,'.xlsx'],input_xlsx_tab_str);
            % input_trajectory_data : [t,theta,r]
            trajectory_t = input_trajectory_data(:,1);
            trajectory_theta = input_trajectory_data(:,2);
            trajectory_r = input_trajectory_data(:,3);

            t_array_rem = rem(t_array, max(trajectory_t));  % max(trajectory_t) is period of the gait
            t_array_floor = floor(t_array/ max(trajectory_t)); % period number

            Fn_traj_theta = @(t) interp1(trajectory_t, (trajectory_theta)/180*pi, t,'linear','extrap');
            theta_array = Fn_traj_theta(t_array_rem);

            % shift theta to the origin to reduce omega spike 
            theta_array_trimed = theta_array -(theta_array(1) - (theta_array(2) - theta_array(1)));
            theta_array_full = max(theta_array_trimed)*t_array_floor + theta_array_trimed;

            theta_array_shifted_ind = find(theta_array_full >= theta_initial_assigned_rad,1,'first');
            theta_array_full_shifted = circshift(theta_array_full, -(theta_array_shifted_ind-1));
            theta_array_full_shifted(end-theta_array_shifted_ind+2:end) = ...
                        theta_array_full_shifted(end-theta_array_shifted_ind+2:end) + theta_array_full(end); 

            Fn_traj_r = @(t) interp1(trajectory_t, trajectory_r, t,'linear','extrap');
            r_array = Fn_traj_r(t_array_rem);
            r_array_shifted = circshift(r_array, -(theta_array_shifted_ind-1));
            trajectory.leg_inertia = 0.02211; %leg
    end

    trajectory.table = [ t_array;
                        theta_array_full_shifted;
                        r_array_shifted];
    trajectory.hip_ini = hip_joint_initial;


    data = Fn_rolling_simulation_sin(landscpae_var , trajectory);
    record_data(record_data_row_index, :) ...
                     = [landscpae_var.amp, ...
                        landscpae_var.freq, ...
                        data.landscape_analysis,...  % [mean, std, base_freq, peak2peak]   % standard deviation is the square root of variance
                        trajectory_mode,...
                        theta_initial_assigned_deg,...
                        data.hip_joint_vs_landscape_length_ratio,...
                        data.work_per_landscape_length,...
                        data.average_speed_x,...
                        data.hip_joint_y_delta];

    record_data_row_index = record_data_row_index + 1;
end  % end of varing initial theta

%% one point data analysis, data of different initial value

analysis_data = [record_data(1,1:2),record_data(1,7),...
                mean(record_data(:,9)),std(record_data(:,9)),...  % hip_joint_vs_landscape_length_ratio
                mean(record_data(:,10)),std(record_data(:,10)),...  % work_per_landscape_length
                mean(record_data(:,11)),std(record_data(:,11)),...  % average_speed_x
                mean(record_data(:,12)),std(record_data(:,12))];  % hip_joint_y_delta_sum







%% 
% subplot(3,1,1)
% plot(record_data(1:20,2),record_data(1:20,9),'linewidth',1.5);
% hold on;
% plot(record_data(21:40,2),record_data(21:40,9),'linewidth',1.5);
% title('freq vs work per length','fontsize', 12);
% xlabel('freq');
% ylabel('work per landscape length');
% legend('dr = 0','dr = 0.045');
% 
% subplot(3,1,2)
% plot(record_data(1:20,2),record_data(1:20,10),'linewidth',1.5);
% hold on;
% plot(record_data(21:40,2),record_data(21:40,10),'linewidth',1.5);
% title('freq vs traveled ratio','fontsize', 12);
% xlabel('freq');
% ylabel('hip/landscape traveled ratio');
% legend('dr = 0','dr = 0.045');
% 
% subplot(3,1,3)
% plot(record_data(1:20,2),record_data(1:20,11),'linewidth',1.5);
% hold on;
% plot(record_data(21:40,2),record_data(21:40,11),'linewidth',1.5);
% title('freq vs average speed','fontsize', 12);
% xlabel('freq');
% ylabel('average speed');
% legend('dr = 0','dr = 0.045');




% 
% figure(1)
% 
% subplot(3,1,1)
% plot(record_data(1:20,2),record_data(1:20,9),'linewidth',1.5);
% hold on;
% plot(record_data(21:40,2),record_data(21:40,9),'linewidth',1.5);
% title('amp vs work per length','fontsize', 12);
% xlabel('amp');
% ylabel('work per landscape length');
% legend('dr = 0','dr = 0.045');
% 
% subplot(3,1,2)
% plot(record_data(1:20,2),record_data(1:20,10),'linewidth',1.5);
% hold on;
% plot(record_data(21:40,2),record_data(21:40,10),'linewidth',1.5);
% title('amp vs traveled ratio','fontsize', 12);
% xlabel('amp');
% ylabel('hip/landscape traveled ratio');
% legend('dr = 0','dr = 0.045');
% 
% subplot(3,1,3)
% plot(record_data(1:20,2),record_data(1:20,11),'linewidth',1.5);
% hold on;
% plot(record_data(21:40,2),record_data(21:40,11),'linewidth',1.5);
% title('amp vs average speed','fontsize', 12);
% xlabel('amp');
% ylabel('average speed');
% legend('dr = 0','dr = 0.045');


% txt = ['T=',num2str(t_end ),'(s), '...
%       'Theta=',num2str(theta_initial*180/pi),'~',num2str(theta_end*180/pi),'(deg), '...
%       'y = A*sin(freq*x)+b, ',...
%       '\mu_s = ', num2str(mu_s),...
%       ', \mu_k = ', num2str(mu_k)]; 
% text( 0.05 , 0.10, txt,'color', 'k','fontsize', 12);
