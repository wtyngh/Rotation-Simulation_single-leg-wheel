% command_type = 'ab';
% 
% 
% switch(command_type)
%     case 'ab'
%         disp("ab");
%     case 'cd'
%         disp("cd");
%     otherwise
%         
% end
        
%%
clear variables; clc;

%% 

% a = linspace(1,1,100);
x = 0:0.1:10;
y = gaussmf(x,[2 5]);
y1 = gaussmf(x,[3 5]);
plot(x,y)
hold on;
plot(x,y1)
xlabel('gaussmf')
% mean(y)
% var(y)


%% Generate random data with given mean and variance
clear variables; clc; close all;
% rng(0,'twister');
a = 5;  %  std = a ; variance = a^2 
b = 500;  % mean = b
c = 3.326;  % frequency
data_y = (a.*randn(1000,1) + b)';


data_x = linspace(0,10,1000);
% data_y = a*sin(c*data_x) + b;

stats = [mean(data_y) std(data_y) var(data_y)]

f = fit(data_x',data_y','fourier2')
f.w
% data_y_sort = sort(data_y);
% 
% stats_sort = [mean(data_y_sort) std(data_y_sort) var(data_y_sort)]


figure
plot(data_y);
% hold on; 
% plot(data_y_sort);
% 
% figure
% histogram(data_y)


%%
a = 2;  % variance = a^2 ; std = a
b = 500;  % mean = b

x_pdf = linspace(-10,10,100);
y_pdf = pdf('Normal',x_pdf, 0, 1); % mu , sigma
plot(x_pdf,y_pdf);

data_y = y_pdf + b;

stats = [mean(data_y) std(data_y) var(data_y)]

% pd = fitdist(data,'Normal')
% pd_test = makedist('Normal',500,5);




% [p,S,mu] = polyfit(data_x,data_y,5);
% poly_data_y = polyval(p,data_x);

% figure
% plot(data_x , data_y);
% hold on;
% plot(data_x , poly_data_y);







%%  write header to excel
% data=ones(10,4);      %Sample 2-dimensional data
% data_cells=num2cell(data);
% col_header={'Temperature','Pressure','X','Y'};     %Row cell array (for column labels)
% row_header(1:10,1)={'Time'};     %Column cell array (for row labels)
% data_out = [col_header;data_cells];
% % data_out_mat = cell2mat( data_out);
% xlswrite('My_file.xls',data,'Sheet1','A2');     %Write data
% xlswrite('My_file.xls',col_header,'Sheet1','B1');     %Write column header
% xlswrite('My_file.xls',row_header,'Sheet1','A2');      %Write row header
% %%
% data=ones(10,4);     %Sample 2-dimensional data
% data_cells=num2cell(data);     %Convert data to cell array
% col_header={'Temperature','Pressure','X','Y'};     %Row cell array (for column labels)
% row_header(1:10,1)={'Time'};     %Column cell array (for row labels)
% output_matrix=[{' '} col_header; row_header data_cells];     %Join cell arrays
% xlswrite('My_file.xls',output_matrix);     %Write data and both headers
% %%
% 
% T = table(['M';'F';'M'],[45 45;41 32;40 34],...
%     {'NY';'CA';'MA'},[true;false;false])
% writetable(T,'myData.xls','Sheet',2,'Range','B2:F6');

