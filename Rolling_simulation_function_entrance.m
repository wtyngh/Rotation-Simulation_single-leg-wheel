%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is used for observing the difference between leg and wheel
% Especially for the different characteristics on different terrains
%
% Geometry included
% Dynamic condisered
%
% Last advised : 2018/02/12
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% opengl info


%% Draw the continuous animation with given conditions 
clear variables; clc;
timer_total = tic;


% 1:torque control, input 'time-torque' trajectory
% 2:Position control, input 'time-theta' trajectory


leg_mass = 1 ; % define the mass of the structure (kg)
leg_inertia = 1; % define the inertia of the structure (kg*m^2)

mu_s = 0.9; % define the equivalent static friction constant between the wheel and the ground 
mu_k = 0.8; % define the equivalent dynamic friction constant


for j = 0:1:1

    %% Inital values 
    hip_joint_initial = [0,0.2];  % initail position of the hip joint
    theta_initial = 0; % define the intial posture of the leg
    theta_end = theta_initial + 4 * pi; % define the fianl posture of the leg
    if j == 0
        delta_r_initial = 0;
    else
        delta_r_initial = 0.045;
    end

    % define how much time the leg is going to run (sec)
    t_initial = 0;  % (s)
    t_end = 10; 

    % define the resolution of the animation
    % More points, higher resolution 
    num_of_iterations = 501;

    t_array = linspace(t_initial, t_end, num_of_iterations);  % t
    % define the gait table
    theta_array = linspace(theta_initial, theta_end, num_of_iterations); % constant omega 
    r_array = 0 * theta_array + delta_r_initial ;  % constant delta_r

    trajectory_table = [ t_array;
                         theta_array;
                         r_array];


    %% Define landscape 

    landscpae_var.x_range = [-0.2, 1.8]; % range of the window
    landscpae_var.y_range = [-0.2, 0.6];
    landscpae_var.x_partition_diff = 0.001; % define the resolution of the gound

    % amp * sin(freq * x) + bias

    for i = 1:1:20

        landscpae_var.amp = 0.01*(i-1);
        landscpae_var.freq = 13;
        landscpae_var.bias = 0;

        try
            data = Fn_rolling_simulation(landscpae_var , trajectory_table);
            record_data(j*20+i,:) = [delta_r_initial,...
                                landscpae_var.amp, landscpae_var.freq, landscpae_var.bias,...
                                data.landscape_analysis,...  % mean, std, var, base_freq
                                data.work_per_landscape_length,...
                                data.hip_joint_vs_landscape_length_ratio...
                                data.average_speed];
        catch
        end
    end
end
% 
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

figure(1)

subplot(3,1,1)
plot(record_data(1:20,2),record_data(1:20,9),'linewidth',1.5);
hold on;
plot(record_data(21:40,2),record_data(21:40,9),'linewidth',1.5);
title('amp vs work per length','fontsize', 12);
xlabel('amp');
ylabel('work per landscape length');
legend('dr = 0','dr = 0.045');

subplot(3,1,2)
plot(record_data(1:20,2),record_data(1:20,10),'linewidth',1.5);
hold on;
plot(record_data(21:40,2),record_data(21:40,10),'linewidth',1.5);
title('amp vs traveled ratio','fontsize', 12);
xlabel('amp');
ylabel('hip/landscape traveled ratio');
legend('dr = 0','dr = 0.045');

subplot(3,1,3)
plot(record_data(1:20,2),record_data(1:20,11),'linewidth',1.5);
hold on;
plot(record_data(21:40,2),record_data(21:40,11),'linewidth',1.5);
title('amp vs average speed','fontsize', 12);
xlabel('amp');
ylabel('average speed');
legend('dr = 0','dr = 0.045');


txt = ['T=',num2str(t_end ),'(s), '...
      'Theta=',num2str(theta_initial*180/pi),'~',num2str(theta_end*180/pi),'(deg), '...
      'y = A*sin(freq*x)+b, ',...
      '\mu_s = ', num2str(mu_s),...
      ', \mu_k = ', num2str(mu_k)]; 
text( 0.05 , 0.10, txt,'color', 'k','fontsize', 12);
