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

%% Inital values

a = [1,1];
b = [2,2];
linspace(1,2)




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

