clear variables;clc;close all;

x_range = [-0.5, 4.5]; % range of the window
y_range = [-0.25, 0.8];

x_partition_diff = 0.001; % define the resolution of the gound
x_partition = x_range(1):x_partition_diff:x_range(2);  % x_partition

%**********************************
freq_var_inc = 0.01;
freq_var_array = 0:freq_var_inc:0.8;
size_freq_array = size(freq_var_array,2);

%***********************************
amp_var_inc = 0.01; 
amp_var_array = 0 : amp_var_inc : 1.0; %[m]
size_amp_array = size(amp_var_array,2);
%*********************************** 
landscape_IO_compare_sin = double.empty(0,3);
landscape_IO_compare_data_row = 1;

freq = 0.2;
amp = 0.2;
bias = 0;


for amp = amp_var_array(1:end)

    landscape_function = @(x) amp * sin( freq *2*pi *x ) + bias ;
    landscape_partition = landscape_function(x_partition);

    % add noise to data
    landscape_partition_noise = awgn(landscape_partition,25,'measured');
    %  20 dB signal-to-noise ratio (SNR) 
    
%     landscape_partition_noise = amp*0.05*rand(size(x_partition)) + landscape_partition;
    
%     plot(x_partition,landscape_partition);
%     hold on;
%     plot(x_partition,landscape_partition_noise);
%     legend('orginal','20 dB SNR');
%     title('Add noise to signal');
%     xlabel('x');
%     ylabel('y');
    

    %=== Retrive landscape freq ===
    % landscape_fourier_fit = fit(x_partition',landscape_partition','fourier2');
    FFT_freq_sample = 20; % Hz % Sampling frequency                           
    FFT_L = 1500;             % Length of signal
    FFT_t = (0:FFT_L-1)*(1/FFT_freq_sample);        % Time vector


    frequecy_response_array = find_dominate_freq(landscape_partition_noise', FFT_freq_sample);

    % mean = mean(landscape_partition)
    std_resp = std(landscape_partition_noise);
    freq_resp = frequecy_response_array(1,1);
    peak2peak_resp = peak2peak(landscape_partition_noise);

%     landscape_IO_compare_sin(landscape_IO_compare_data_row,:) = [freq,freq_resp];
    landscape_IO_compare_sin(landscape_IO_compare_data_row,:) = [amp,std_resp,peak2peak_resp];

    landscape_IO_compare_data_row = landscape_IO_compare_data_row + 1;
end


%%
[curve_fitted, fit_info] = fit(landscape_IO_compare_sin(:,2),landscape_IO_compare_sin(:,1),'poly2');

%%
plot(landscape_IO_compare_sin(:,2),landscape_IO_compare_sin(:,1));
title('Data std vs amp input');
xlabel('Data std');
ylabel('Amp input');

plot_x_range = xlim;
plot_y_range = ylim;

text(plot_x_range(1)+0.1,plot_y_range(2)-0.1,'(Real freq) = 1.412*(Freq analysis)','fontsize',12);
text(plot_x_range(1)+0.1,plot_y_range(2)-0.15, 'r^2 = 0.9999','fontsize',12);



