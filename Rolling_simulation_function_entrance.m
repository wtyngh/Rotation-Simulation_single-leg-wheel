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




%% Inital values 
hip_joint_initial = [0,0.2];  % initail position of the hip joint
theta_initial = 0; % define the intial posture of the leg
theta_end = theta_initial + 2 * pi; % define the fianl posture of the leg
delta_r_initial = 0.045;

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

landscpae_var.x_range = [-0.2, 1.5]; % range of the window
landscpae_var.y_range = [-0.2, 0.6];
landscpae_var.x_partition_diff = 0.001; % define the resolution of the gound

% amp * sin(freq * x) + bias
landscpae_var.amp = 0.1;
landscpae_var.freq = 10;
landscpae_var.bias = 0;

data = Fn_rolling_simulation(landscpae_var , trajectory_table);

