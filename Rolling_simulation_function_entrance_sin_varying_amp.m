%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This code is used for observing the difference between leg and wheel
% Especially for the different characteristics on different terrains
% For exploring data points
%
% Varying amplitude
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


leg_mass = 5 ; % define the mass of the structure (kg)
% define the inertia of the structure (kg*m^2) predict by solidworks
% wheel:0.1118(kg*m^2)
% leg:0.02211(kg*m^2)
% motor inertia is small enough to be neglected

mu_s = 0.6; % define the equivalent static friction constant between the wheel and the ground 
mu_k = 0.5; % define the equivalent dynamic friction constant

record_data_row_index = 1;
analysis_data_row_index = 1; % index for recording one point, with varying initial theta
%%
enable.savedata = 1; 
enable.savefig = 1;

%% define how much time the leg is going to run (sec)
t_initial = 0;  % (s)
t_end = 5; 
% define the resolution of the animation
% More points, higher resolution 
num_of_iterations = 1001;
t_array = linspace(t_initial, t_end, num_of_iterations);  % t

% initialize data_record
record_data = double.empty(0,14);


%% Define landscape 

landscpae_var.x_range = [-1.2, 4.5]; % range of the window
landscpae_var.y_range = [-0.25, 0.8];
landscpae_var.x_partition_diff = 0.001; % define the resolution of the gound

% amp * sin(freq * x) + bias

% for amp_var_ind = 1:1:amp_var_size
% for amp_var_ind = 1:1:amp_var_size

%**********************************
landscpae_var.mu_s = mu_s;
landscpae_var.mu_k = mu_k;

%**********************************
freq_var_inc = 0.1;
freq_var_array = 0.9:freq_var_inc:1.0;
size_freq_array = size(freq_var_array,2);

%***********************************
amp_var_inc = 0.02; 
amp_var_array = 0.02 : amp_var_inc : 0.2; %[m]
size_amp_array = size(amp_var_array,2);
%***********************************

%% Inital values 

hip_joint_initial = [0,0.2];  % initail position of the hip joint    


theta_initial_assigned_deg_array = 0:10:179;  % 18 points
% theta_initial_assigned_deg_array = 10;
size_theta_initial_array = size(theta_initial_assigned_deg_array,2);

total_loop_iteration = size_freq_array * size_amp_array * size_theta_initial_array * 4;
% initialize analysis_data
analysis_data = double.empty(0,13);

%% import xlsx data

input_trajectory_data_filename = 'CPG trajectory';
trot_xlsx_tab_str = 'Trot, V=400';
trot_trajectory_data = xlsread([input_trajectory_data_filename,'.xlsx'],trot_xlsx_tab_str);
walk_xlsx_tab_str = 'Walk, V=400';
walk_trajectory_data = xlsread([input_trajectory_data_filename,'.xlsx'],walk_xlsx_tab_str);

%% Start main loop
timer_total = tic;

for freqency = freq_var_array(1:end)

    for trajectory_mode = 1:4  
    % 1:constant omega, constant dr = 0, wheel mode 
    % 2:constant omega, constant dr = 0.045, legged mode 
    % 3:assigned trajectory, CPG trajectory is used here, Trot
    % 4:assigned trajectory, CPG trajectory is used here, Walk
        for amplitude = amp_var_array(1:end)
            landscpae_var.amp = amplitude;
            landscpae_var.freq = freqency;
            landscpae_var.bias = 0;


            for theta_initial_assigned_deg = theta_initial_assigned_deg_array(1:end)

                theta_initial_assigned_rad = theta_initial_assigned_deg/180*pi;

                timer_loop = tic;

                % define the gait table                    
                switch trajectory_mode

                    case {1,2} % constant omega, constant dr
                        forward_vel_set = 0.4; % define forward velocity
                        if trajectory_mode == 1 % wheel mode
%                             hip_joint_initial = [0,0.1];  % bug fixing for wheel critical condition
                            delta_r_initial = 0;
                            forward_dis = forward_vel_set / 0.11*(t_end-t_initial);
                            theta_end = theta_initial_assigned_rad + forward_dis; %(V/r)*t=w*t
                            trajectory.leg_inertia = 0.1118; % wheel
                        else  % trajectory_mode == 2 % legged mode
                            hip_joint_initial = [0,0.2];  % bug fixing for wheel critical condition
                            delta_r_initial = 0.045;
                            forward_dis = forward_vel_set / 0.155*(t_end-t_initial) * 1.2; % fixing constant
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
%                         hip_joint_initial = [0,0.2];  % bug fixing for wheel critical condition
                        % Load trajectory data                    
                        input_trajectory_data_filename = 'CPG trajectory';
                        if trajectory_mode == 3
                            input_xlsx_tab_str = trot_xlsx_tab_str;
                            input_trajectory_data = trot_trajectory_data;
                        else  % trajectory_mode == 4
                            input_xlsx_tab_str = walk_xlsx_tab_str;
                            input_trajectory_data = walk_trajectory_data;
                        end
                        trajectory.name = [input_trajectory_data_filename,', ',input_xlsx_tab_str,'[mm/s]'];
%                         input_trajectory_data = xlsread([input_trajectory_data_filename,'.xlsx'],input_xlsx_tab_str);

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
                                    mu_s,mu_k,...
                                    data.landscape_analysis,...  % [mean, std, base_freq, peak2peak]   % standard deviation is the square root of variance
                                    trajectory_mode,...
                                    theta_initial_assigned_deg,...
                                    data.hip_joint_vs_landscape_length_ratio,...
                                    data.work_per_landscape_length,...
                                    data.average_speed_x,...
                                    data.hip_joint_y_delta];

                time_str = [sprintf('%.1f',(record_data_row_index/total_loop_iteration*100)),'%% , ',...
                            sprintf('loop = %.2f(s)', toc(timer_loop)),' , ',...
                            sprintf('Elapsed = %.2f(s)\n', toc(timer_total))];
        %                     sprintf('%.1f',(analysis_data_row_index/size_amp_array*100)),'%% , ',...

                fprintf(time_str);                


                record_data_row_index = record_data_row_index + 1;
            end  % end of varing initial theta

            %% one point data analysis, data of different initial value

            record_data_row_range = ... % different initial theta points are combined as a set
                (analysis_data_row_index-1) * size_theta_initial_array + 1: analysis_data_row_index * size_theta_initial_array;

            analysis_data(analysis_data_row_index,:) = ...
                           [record_data(record_data_row_range(1),1:4),record_data(record_data_row_range(1),9),... % [amp, freq, mu_s, mu_k, mode]
                            mean(record_data(record_data_row_range,11)),std(record_data(record_data_row_range,11)),...  % hip_joint_vs_landscape_length_ratio
                            mean(record_data(record_data_row_range,12)),std(record_data(record_data_row_range,12)),...  % work_per_landscape_length
                            mean(record_data(record_data_row_range,13)),std(record_data(record_data_row_range,13)),...  % average_speed_x
                            mean(record_data(record_data_row_range,14)),std(record_data(record_data_row_range,14))];  % hip_joint_y_delta_sum

            analysis_data_row_index = analysis_data_row_index + 1;

        end  % end of varying amplitude

    end  % end of varying trajectory mode

    %% Save analysis data, for fixed freq, varing amp

    if enable.savedata == 1
        data_file_name = ['Fs=',num2str(mu_s),',Fk=',num2str(mu_k),',Amp[',num2str(amp_var_array(1)),'~',num2str(amp_var_array(end)),...
                    '],freq=',num2str(freqency),'.mat'];
        save_file_path = fullfile(pwd,'Sim result_sin_400_adj_fric',data_file_name);

        save(save_file_path,'analysis_data');
    end

    
    
    



%%
% %% Crossing ratio plot
% 
% figure(2)
% hold on;
% for mode = 1:4  
%     selected_data_row_ind = find(analysis_data(:,5) == mode);
%     % close all;
%     errorbar(analysis_data(selected_data_row_ind,1),analysis_data(selected_data_row_ind,6),analysis_data(selected_data_row_ind,7));  % crossing ratio
% end
% hold off;
% 
% title('Amplitude vs Corssing ratio, with dif. theta I.C.','fontsize', 12);
% xlabel('Amplitude [m]');
% ylabel('Crossing ratio [m/m]');
% 
% xl = xlim;
% yl = ylim;
% 
% legend('Fixed dr = 0 [m]','Fixed dr = 0.045 [m]','Trot, V=400 [m/s]','Walk, V=400 [m/s]');
% 
% freq_txt = ['Freq = ',num2str(analysis_data(1,2),'%.2f'),' [Hz]'] ;
% mu_s_txt = ['\mu_s = ',num2str(analysis_data(1,3),'%.1f')] ;
% mu_k_txt = ['\mu_k = ',num2str(analysis_data(1,4),'%.1f')] ;
% def_txt = 'Crossing ratio = traj hip / traj landscape';
% 
% 
% text( (xl(2)-xl(1))*0.02 , yl(2)-(yl(2) - yl(1))*0.05 , freq_txt ,'color', 'k', 'fontsize', 10);
% text( (xl(2)-xl(1))*0.02 , yl(2)-(yl(2) - yl(1))*0.1 , mu_s_txt ,'color', 'k', 'fontsize', 10);
% text( (xl(2)-xl(1))*0.02 , yl(2)-(yl(2) - yl(1))*0.15 , mu_k_txt ,'color', 'k', 'fontsize', 10);
% text( (xl(2)-xl(1))*0.02 , yl(2)-(yl(2) - yl(1))*0.2 , def_txt ,'color', 'k', 'fontsize', 10);
% 
% 
% 
% 
% 
% 
% 
% %% Work plot 
% 
% figure(3)
% hold on;
% for mode = 1:4  
%     selected_data_row_ind = find(analysis_data(:,5) == mode);
%     % close all;
%     errorbar(analysis_data(selected_data_row_ind,1),analysis_data(selected_data_row_ind,8),analysis_data(selected_data_row_ind,9));  % work
% end
% hold off;
% 
% title('Amplitude vs Work per length, with dif. theta I.C.','fontsize', 12);
% xlabel('Amplitude [m]');
% ylabel('work per landscape length [J/m]');
% 
% xl = xlim;
% yl = ylim;
% 
% legend('Fixed dr = 0 [m]','Fixed dr = 0.045 [m]','Trot, V=400 [m/s]','Walk, V=400 [m/s]');
% 
% freq_txt = ['Freq = ',num2str(analysis_data(1,2),'%.2f'),' [Hz]'] ;
% mu_s_txt = ['\mu_s = ',num2str(analysis_data(1,3),'%.1f')] ;
% mu_k_txt = ['\mu_k = ',num2str(analysis_data(1,4),'%.1f')] ;
% 
% text( (xl(2)-xl(1))*0.02 , yl(2)-(yl(2) - yl(1))*0.05 , freq_txt ,'color', 'k', 'fontsize', 10);
% text( (xl(2)-xl(1))*0.02 , yl(2)-(yl(2) - yl(1))*0.1 , mu_s_txt ,'color', 'k', 'fontsize', 10);
% text( (xl(2)-xl(1))*0.02 , yl(2)-(yl(2) - yl(1))*0.15 , mu_k_txt ,'color', 'k', 'fontsize', 10);
% 
% 
% %% Vel plot 
% 
% figure(4)
% hold on;
% for mode = 1:4  
%     selected_data_row_ind = find(analysis_data(:,5) == mode);
%     % close all;
%     errorbar(analysis_data(selected_data_row_ind,1),analysis_data(selected_data_row_ind,10),analysis_data(selected_data_row_ind,11));  % work
% end
% hold off;
% 
% title('Amplitude vs avg Vx, with dif. theta I.C.','fontsize', 12);
% xlabel('Amplitude [m]');
% ylabel('average Vx  [m/s]');
% 
% xl = xlim;
% yl = ylim;
% 
% legend('Fixed dr = 0 [m]','Fixed dr = 0.045 [m]','Trot, V=400 [m/s]','Walk, V=400 [m/s]');
% 
% freq_txt = ['Freq = ',num2str(analysis_data(1,2),'%.2f'),' [Hz]'] ;
% mu_s_txt = ['\mu_s = ',num2str(analysis_data(1,3),'%.1f')] ;
% mu_k_txt = ['\mu_k = ',num2str(analysis_data(1,4),'%.1f')] ;
% 
% 
% text( (xl(2)-xl(1))*0.02 , (yl(2) - yl(1))*0.09 , freq_txt ,'color', 'k', 'fontsize', 10);
% text( (xl(2)-xl(1))*0.02 , (yl(2) - yl(1))*0.05 , mu_s_txt ,'color', 'k', 'fontsize', 10);
% text( (xl(2)-xl(1))*0.02 , (yl(2) - yl(1))*0.01 , mu_k_txt ,'color', 'k', 'fontsize', 10);
% 
% %% Height variance plot 
% 
% figure(5)
% hold on;
% for mode = 1:4  
%     selected_data_row_ind = find(analysis_data(:,5) == mode);
%     % close all;
%     errorbar(analysis_data(selected_data_row_ind,1),analysis_data(selected_data_row_ind,12),analysis_data(selected_data_row_ind,13));  % work
% end
% hold off;
% 
% title('Amplitude vs Hip height delta sum, with dif. theta I.C.','fontsize', 12);
% xlabel('Amplitude [m]');
% ylabel('Hip height delta sum  [m]');
% 
% xl = xlim;
% yl = ylim;
% 
% legend('Fixed dr = 0 [m]','Fixed dr = 0.045 [m]','Trot, V=400 [m/s]','Walk, V=400 [m/s]');
% 
% freq_txt = ['Freq = ',num2str(analysis_data(1,2),'%.2f'),' [Hz]'] ;
% mu_s_txt = ['\mu_s = ',num2str(analysis_data(1,3),'%.1f')] ;
% mu_k_txt = ['\mu_k = ',num2str(analysis_data(1,4),'%.1f')] ;
% 
% 
% text( (xl(2)-xl(1))*0.02 , yl(2)-(yl(2) - yl(1))*0.08 , freq_txt ,'color', 'k', 'fontsize', 10);
% text( (xl(2)-xl(1))*0.02 , yl(2)-(yl(2) - yl(1))*0.16 , mu_s_txt ,'color', 'k', 'fontsize', 10);
% text( (xl(2)-xl(1))*0.02 , yl(2)-(yl(2) - yl(1))*0.24 , mu_k_txt ,'color', 'k', 'fontsize', 10);

%% Plot together

    set(gcf,'name','Plot sim result','Position', [100 100 1500 800]);
    subplot(2,2,1);
    hold on;
    for mode = 1:4  
        selected_data_row_ind = find(analysis_data(:,5) == mode);
        % close all;
        errorbar(analysis_data(selected_data_row_ind,1),analysis_data(selected_data_row_ind,6),analysis_data(selected_data_row_ind,7));  % crossing ratio
    end
    hold off;

    title('Amplitude vs Corssing ratio, with dif. theta I.C.','fontsize', 12);
    xlabel('Amplitude [m]');
    ylabel('Crossing ratio [m/m]');

    xl = xlim;
    yl = ylim;

    % legend('Fixed dr = 0 [m]','Fixed dr = 0.045 [m]','Trot, V=400 [m/s]','Walk, V=400 [m/s]');

    % freq_txt = ['Freq = ',num2str(analysis_data(1,2),'%.2f'),' [Hz]'] ;
    % mu_s_txt = ['\mu_s = ',num2str(analysis_data(1,3),'%.1f')] ;
    % mu_k_txt = ['\mu_k = ',num2str(analysis_data(1,4),'%.1f')] ;
    def_txt = 'Crossing ratio = traj hip / traj landscape';


    % text( (xl(2)-xl(1))*0.02 , yl(2)-(yl(2) - yl(1))*0.08 , freq_txt ,'color', 'k', 'fontsize', 10);
    % text( (xl(2)-xl(1))*0.02 , yl(2)-(yl(2) - yl(1))*0.16 , mu_s_txt ,'color', 'k', 'fontsize', 10);
    % text( (xl(2)-xl(1))*0.02 , yl(2)-(yl(2) - yl(1))*0.24 , mu_k_txt ,'color', 'k', 'fontsize', 10);
    text( (xl(2)-xl(1))*0.01 , yl(2)-(yl(2) - yl(1))*0.05 , def_txt ,'color', 'k', 'fontsize', 10); 

    subplot(2,2,2);
    hold on;
    for mode = 1:4  
        selected_data_row_ind = find(analysis_data(:,5) == mode);
        % close all;
        errorbar(analysis_data(selected_data_row_ind,1),analysis_data(selected_data_row_ind,8),analysis_data(selected_data_row_ind,9));  % work
    end
    hold off;

    title('Amplitude vs Work per length, with dif. theta I.C.','fontsize', 12);
    xlabel('Amplitude [m]');
    ylabel('work per landscape length [J/m]');

    xl = xlim;
    yl = ylim;

    % legend('Fixed dr = 0 [m]','Fixed dr = 0.045 [m]','Trot, V=400 [m/s]','Walk, V=400 [m/s]');
    % 
    % freq_txt = ['Freq = ',num2str(analysis_data(1,2),'%.2f'),' [Hz]'] ;
    % mu_s_txt = ['\mu_s = ',num2str(analysis_data(1,3),'%.1f')] ;
    % mu_k_txt = ['\mu_k = ',num2str(analysis_data(1,4),'%.1f')] ;
    def_txt = 'Work per length = abs work / traveled landscape x';
    % 
    % text( (xl(2)-xl(1))*0.02 , yl(2)-(yl(2) - yl(1))*0.08 , freq_txt ,'color', 'k', 'fontsize', 10);
    % text( (xl(2)-xl(1))*0.02 , yl(2)-(yl(2) - yl(1))*0.16 , mu_s_txt ,'color', 'k', 'fontsize', 10);
    % text( (xl(2)-xl(1))*0.02 , yl(2)-(yl(2) - yl(1))*0.24 , mu_k_txt ,'color', 'k', 'fontsize', 10);
    text( (xl(2)-xl(1))*0.01 , yl(2)-(yl(2) - yl(1))*0.05 , def_txt ,'color', 'k', 'fontsize', 10); 

    subplot(2,2,3);
    hold on;
    for mode = 1:4  
        selected_data_row_ind = find(analysis_data(:,5) == mode);
        % close all;
        errorbar(analysis_data(selected_data_row_ind,1),analysis_data(selected_data_row_ind,10),analysis_data(selected_data_row_ind,11));  % work
    end
    hold off;

    title('Amplitude vs avg Vx, with dif. theta I.C.','fontsize', 12);
    xlabel('Amplitude [m]');
    ylabel('average Vx  [m/s]');

    % xl = xlim;
    % yl = ylim;

    % legend('Fixed dr = 0 [m]','Fixed dr = 0.045 [m]','Trot, V=400 [m/s]','Walk, V=400 [m/s]');
    % 
    % freq_txt = ['Freq = ',num2str(analysis_data(1,2),'%.2f'),' [Hz]'] ;
    % mu_s_txt = ['\mu_s = ',num2str(analysis_data(1,3),'%.1f')] ;
    % mu_k_txt = ['\mu_k = ',num2str(analysis_data(1,4),'%.1f')] ;
    % 
    % text( (xl(2)-xl(1))*0.02 , (yl(2) - yl(1))*0.24 , freq_txt ,'color', 'k', 'fontsize', 10);
    % text( (xl(2)-xl(1))*0.02 , (yl(2) - yl(1))*0.16 , mu_s_txt ,'color', 'k', 'fontsize', 10);
    % text( (xl(2)-xl(1))*0.02 , (yl(2) - yl(1))*0.08 , mu_k_txt ,'color', 'k', 'fontsize', 10);

    subplot(2,2,4);
    hold on;
    for mode = 1:4  
        selected_data_row_ind = find(analysis_data(:,5) == mode);
        % close all;
        errorbar(analysis_data(selected_data_row_ind,1),analysis_data(selected_data_row_ind,12),analysis_data(selected_data_row_ind,13));  % work
    end
    hold off;

    title('Amplitude vs Hip height delta per length, with dif. theta I.C.','fontsize', 12);
    xlabel('Amplitude [m]');
    ylabel('Hip height delta per length [m/m]');

    xl = xlim;
    yl = ylim;

    legend('Fixed dr = 0 [m]','Fixed dr = 0.045 [m]','Trot, V=0.4 [m/s]','Walk, V=0.4 [m/s]');

    freq_txt = ['Freq = ',num2str(analysis_data(end,2),'%.2f'),' [Hz]'] ;
    mu_s_txt = ['\mu_s = ',num2str(analysis_data(end,3),'%.1f')] ;
    mu_k_txt = ['\mu_k = ',num2str(analysis_data(end,4),'%.1f')] ;
    def_txt = 'Hip height delta per length = Hip height delta sum / traveled landscape x';

    text( xl(2)-(xl(2)-xl(1))*0.2 , (yl(2) - yl(1))*0.08 , freq_txt ,'color', 'k', 'fontsize', 10);
    text( xl(2)-(xl(2)-xl(1))*0.2 , (yl(2) - yl(1))*0.24 , mu_s_txt ,'color', 'k', 'fontsize', 10);
    text( xl(2)-(xl(2)-xl(1))*0.2 , (yl(2) - yl(1))*0.16 , mu_k_txt ,'color', 'k', 'fontsize', 10);
    text( (xl(2)-xl(1))*0.01 , yl(2)-(yl(2) - yl(1))*0.05 , def_txt ,'color', 'k', 'fontsize', 10); 

    %% save fig
    if enable.savefig == 1
    
        fig_file_name = ['Fs=',num2str(mu_s),',Fk=',num2str(mu_k),',Amp[',num2str(amp_var_array(1)),'~',num2str(amp_var_array(end)),...
                    '],freq=',num2str(freqency),'.fig'];

        save_fig_path = fullfile(pwd,'Sim result_sin_400_adj_fric',fig_file_name);        

        saveas(gca, save_fig_path);


        close all; 
        clear analysis_data;
    
    end

end

