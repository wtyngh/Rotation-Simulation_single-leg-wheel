function frequecy_response_array = find_dominate_freq(path_points, FFT_Fs)

length_of_path_points =  size(path_points,1);
fft_n = 2^nextpow2(length_of_path_points); % to improve fft performance

FFT_result = fft(path_points, fft_n);

% Calculate the double-sided spectrum and single-sided spectrum of each signal.
P2 = abs(FFT_result/fft_n);
P1 = P2(1:fft_n/2+1, 1);
P1(2:end-1) = 2*P1(2:end-1);

response_axis = P1;
frequency_axis = FFT_Fs*(0:(fft_n/2))/fft_n;


% plot(frequency_axis,response_axis);

adjacency_feq_bond = 1;
cutof_response =  mean(response_axis) + std(response_axis); % user def


sub_sample_grid_num = FFT_Fs/2/10;

% frequecy_response_array = double.empty(6,0);

frequecy_response_array_ind = 1;

for sample_region_ind = 1:sub_sample_grid_num
    freq_low_bound = (sample_region_ind-1) * (FFT_Fs/2)/sub_sample_grid_num;
    freq_high_bound = (sample_region_ind) * (FFT_Fs/2)/sub_sample_grid_num ;
    
    sample_ind =( frequency_axis >= freq_low_bound ) &  (frequency_axis < freq_high_bound);
    [M,I] = max(response_axis(sample_ind));

    freq_res = frequency_axis(I) + freq_low_bound;
    
    if sample_region_ind > 1 
            if abs(frequecy_response_array(frequecy_response_array_ind-1,1) - freq_res) <= adjacency_feq_bond ...
                   % choose the larger response 
                    if M > frequecy_response_array(frequecy_response_array_ind-1,2)
                        frequecy_response_array(frequecy_response_array_ind-1,:) = [freq_res, M];
                        frequecy_response_array_ind = frequecy_response_array_ind + 1;
                    else
                        % do nothing
                    end
            else % non-adacency freq
                if M > cutof_response

                    frequecy_response_array(frequecy_response_array_ind,:) = [freq_res, M];
                    frequecy_response_array_ind = frequecy_response_array_ind + 1;
                end
            end
    else
        frequecy_response_array(sample_region_ind,:) = [freq_res, M];
        frequecy_response_array_ind = frequecy_response_array_ind + 1;
    end
end

% sortrows(frequecy_response_array,2,'descend') % Sort the rows based on the values

end