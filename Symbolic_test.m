% function Symbolic_test (initial_values)


clear variables; clc;
%% define trajectory

radius = 0.11;

% define how much time the leg is going to run (sec)
t_initial = 0;
t_end = 10; 

theta_initial = pi/4; % define the intial posture of the leg
theta_end = theta_initial + 2 * pi; % define the fianl posture of the leg


% define the resolution of the animation
% More points, higher resolution 
num_of_iterations = 1001;


t_array = linspace(t_initial, t_end, num_of_iterations); 
% theta_array = linspace(theta_initial, theta_end, num_of_iterations);
theta_array = linspace(theta_initial, theta_end, num_of_iterations);
r_array = 0 * t_array   ; 


t_increment = (t_end - t_initial)/ (num_of_iterations - 1);
E = 10^-8;

% % for testing
% trajectory_table = [ t_array;
%                      theta_array;
%                      r_array ];


%% Define landscape
x_range = [-0.2, 1.5]; % range of the window
y_range = [-0.2, 0.6];

x_partition_diff = 0.001; % define the resolution of the gound
x_partition = x_range(1):x_partition_diff:x_range(2);  % x_partition

% landscape_partition =  x_partition .^ 2  ;
landscape_partition =  x_partition .* 0 + x_partition ;

landscape_partition_diff = [0, diff(landscape_partition)];
landscape_partition_slope = atan2(landscape_partition_diff, x_partition_diff);

land_table = [ x_partition;
               landscape_partition_diff;
               landscape_partition_slope ];

%% define trajectory function
% syms x t

% landscape function
phi = @(x) interp1(x_partition, landscape_partition_slope, x);
% theta trajectory function
theta = @(t) interp1(t_array,theta_array,t);
% r trajectory function
dr = @(t) interp1(t_array,r_array,t);

%% define geometry

% delta theta, initial position at 3/2 pi
th = @(t) [0,0,-theta(t)]; 

% from half cercle center to contact point
% reference x is with respect to the contact point
R = @(x) radius*[sin(phi(x)),-cos(phi(x)),0]; 

% from hip to the center half circle_1
r1 = @(t) dr(t)*[sin(theta(t)),cos(theta(t)),0];
% from hip to the center half circle_2
r2 = @(t) -dr(t)*[sin(theta(t)),cos(theta(t)),0];  


% Position_1 = @(x,t) cross(R(x),th(t)) - ( R(x) + r1(t) );
Position_1 = @(x,t)  - ( R(x) + r1(t) );
Vel_1 = @(x,t) (Position_1(x,t+t_increment)-Position_1(x,t))/t_increment ...
        + (Position_1(x+x_partition_diff,t)-Position_1(x,t))/x_partition_diff ...
        + cross(R(x), (th(t+t_increment)-th(t))/t_increment ); 
    
    
Acc_1 = @(x,t) (Vel_1(x,t+t_increment)-Vel_1(x,t))/t_increment ...
      + (Vel_1(x+x_partition_diff,t)-Vel_1(x,t))/x_partition_diff;
%%


Position_2 = @(x,t) - ( R(x) + r2(t) );
% velocity: differentail of position
% d/dx + d/dt
Vel_2 = @(x,t) (Position_2(x,t+t_increment)-Position_2(x,t))/t_increment ...
        + (Position_2(x+x_partition_diff,t)-Position_2(x,t))/x_partition_diff ...
        + cross(R(x), (th(t+t_increment)-th(t))/t_increment ); 
    
    
Acc_2 = @(x,t) (Vel_2(x,t+t_increment)-Vel_2(x,t))/t_increment ...
      + (Vel_2(x+x_partition_diff,t)-Vel_2(x,t))/x_partition_diff;


% matlabFunction(simplify(Acc_2), 'file', 'define_Acceraltion', 'vars', {x,t});
% functions.leg_1.Position = Position_1;
% functions.leg_1.Vel = Vel_1;
% functions.leg_1.Acc = Acc_1;
% 
% functions.leg_2.Position = Position_2;
% functions.leg_2.Vel = Vel_2;
% functions.leg_2.Acc = Acc_2;

