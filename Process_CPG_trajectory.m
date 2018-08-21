%%

input_data_filename = '20180409_Wheel, V=200_1';


output_data_filename = 'CPG_trajectory';
xlsx_tab_str = 'Wheel, V=200';

%%
input_data = xlsread([input_data_filename,'.csv']);


%%
trimed_data = input_data(input_data(:,1) ~= 0,:);

t_array = trimed_data(:,1) - trimed_data(1,1);
theta_array = trimed_data(:,2);
R_array = trimed_data(:,3);

%%
subplot(2,1,1)
plot(t_array,theta_array);
subplot(2,1,2)
plot(t_array,R_array);



%%
output_data = [t_array,theta_array,R_array];
output_data(end,:) = [];


[xls_status, xls_message] = xlswrite([output_data_filename,'.xlsx'],output_data, xlsx_tab_str);
%     [xls_status, xls_message] = writetable(data_table,'table.xlsx','Sheet', xlsx_tab_str);
if xls_status == 1
    fprintf('xlsx write sucessful\n');
else
    fprintf('xlsx write error\n');
end