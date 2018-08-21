clear variables; close all; clc;

input_data_filename = 'CPG trajectory';

num_of_iterations = 1000;

for mode = 1:4
    switch mode 
        case 1
            forward_vel_set = 0.4; % define forward velocity
            delta_r_initial = 0;
            theta_array = linspace(0, 360, num_of_iterations); % constant omega
            omega = forward_vel_set / 0.11 ; % fixing constant %(V/r)=w
            t_elapsed = (2*pi)/omega; 
            t_array = linspace(0,t_elapsed,num_of_iterations);
            r_array = 0 * theta_array + delta_r_initial ;  % constant delta_r
        case 2
            forward_vel_set = 0.4; % define forward velocity
            delta_r_initial = 0.045;
            theta_array = linspace(0, 360, num_of_iterations); % constant omega
            omega = forward_vel_set / 0.155 * 1.2; % fixing constant %(V/r)=w
            t_elapsed = (2*pi)/omega; 
            t_array = linspace(0,t_elapsed,num_of_iterations);
            r_array = 0 * theta_array + delta_r_initial ;  % constant delta_r
        case 3
            xlsx_tab_str = 'Trot, V=400';
            input_data = xlsread([input_data_filename,'.xlsx'],xlsx_tab_str);

            t_array = input_data(:,1);
            theta_array = input_data(:,2);
            r_array = input_data(:,3);
        case 4
            xlsx_tab_str = 'Walk, V=400';
            input_data = xlsread([input_data_filename,'.xlsx'],xlsx_tab_str);

            t_array = input_data(:,1);
            theta_array = input_data(:,2);
            r_array = input_data(:,3);
    end


    
    subplot(2,1,1)
    plot(t_array,theta_array, 'linewidth', 2);
    hold on;
    
    subplot(2,1,2)
    plot(t_array,r_array, 'linewidth', 2);
    hold on;
end
subplot(2,1,1)
title('Theta, V=400 [mm/s]');
xlabel('t (s)');
ylabel('Theta (deg)');
legend('Wheel','dr=0.045[m]','Trot','Walk');

subplot(2,1,2)
title('R, V=400 [mm/s]');
xlabel('t (s)');
ylabel('R (m)');
legend('Wheel','dr=0.045[m]','Trot','Walk');