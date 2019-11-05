clear variables;clc;close all;
%% plot stair critical point
clear variables;clc;close all;

data = [
0.2	1	0.9	0.090 	0.120 ;
0.2	0.9	0.8	0.085 	0.120 ;
0.2	0.8	0.7	0.085 	0.115 ;
0.2	0.7	0.6	0.085 	0.110 ;
0.2	0.6	0.5	0.085 	0.105 ;

0.4	1	0.9	0.095 	0.130 ;
0.4	0.9	0.8	0.090 	0.130 ;
0.4	0.8	0.7	0.085 	0.125 ;
0.4	0.7	0.6	0.085 	0.120 ;
0.4	0.6	0.5	0.085 	0.115 ;

];

data(:,6) = data(:,5) - data(:,4); %for ploting area 

drawing_mode = 1; 
% 1:draw the two critical point plots of differernt V together 
% 2:draw the critical point plot V = 0.2 m/s and the corresponding point of the stairs 
  
x = data(1:5,2);
x1 = [ fliplr(x') , x'];
x2 = [x', fliplr(x')];

title_fontsize = 7;
axis_fontsize = 7;
legend_fontsize = 6.5;

if drawing_mode == 1
    
    outputsize = [9 5.5]; % centimeters

    subplot(1,2,1)

    plot(x,data(1:5,4),'o-.','linewidth',1.2);
    hold on;
    plot(x,data(1:5,5),'*-','linewidth',1.2);

    Between1 = [zeros(1,size(x,1)) ,data(1:5,4)'];
    f1 = fill(x1, Between1, [ 0.9290    0.6940    0.1250],'EdgeColor','none');
    % set(gca,'Layer','bottom');
    Between2 = [data(1:5,4)', fliplr(data(1:5,5)')];
    f2 = fill(x2, Between2, [0.3010    0.7450    0.9330],'EdgeColor','none');
    % set(gca,'Layer','bottom');
    % set(h,'edgecolor','white');
    h1 = get(gca,'Children');
    set(gca,'Children',[h1(3) h1(4) h1(1) h1(2)]);

    grid on;
    box on;
    xlim ([0.6 1.0]);
    ylim ([0 0.14]);



    title('\mu_s vs L_h , V = 0.2 [m/s]','fontsize',title_fontsize);
    xlabel('\mu_s' ,'fontsize',axis_fontsize);
    ylabel('L_h [m]','fontsize',axis_fontsize);

    set(gca,'Xtick',0.6:0.1:1.0,'Ytick',0:0.02:0.18,'fontsize',axis_fontsize);

    %%
    subplot(1,2,2)

    plot(x,data(6:10,4),'o-.','linewidth',1.2);
    hold on;
    plot(x,data(6:10,5),'*-','linewidth',1.2);

    Between3 = [zeros(1,size(x,1)) ,data(6:10,4)'];
    fill(x1, Between3, [ 0.9290    0.6940    0.1250],'EdgeColor','none');
    % set(gca,'Layer','bottom');
    Between4 = [data(6:10,4)', fliplr(data(6:10,5)')];
    fill(x2, Between4, [0.3010    0.7450    0.9330],'EdgeColor','none');
    % set(h,'edgecolor','white');
    h2 = get(gca,'Children');
    set(gca,'Children',[h2(3) h2(4) h2(1) h2(2)]);

    grid on;
    box on;
    xlim ([0.6 1.0]);
    ylim ([0 0.14]);

    title('\mu_s vs L_h , V = 0.4 [m/s]','fontsize',title_fontsize);
    xlabel('\mu_s' ,'fontsize',axis_fontsize);
    ylabel('L_h [m]','fontsize',axis_fontsize);

    set(gca,'Xtick',0.6:0.1:1.0,'Ytick',0:0.02:0.18,'fontsize',axis_fontsize);

    %%
    lgd = legend({'Wheel OP zone','leg OP zone','L_{h,max,wheeled}','L_{h,max,legged}'},...
                'location','southeast');
    lgd.FontSize = legend_fontsize;


else
    
    outputsize = [4.8 5.2]; % centimeters

    plot(x,data(1:5,4),'o-.','linewidth',1.2);
    hold on;
    plot(x,data(1:5,5),'*-','linewidth',1.2);

    Between1 = [zeros(1,size(x,1)) ,data(1:5,4)'];
    f1 = fill(x1, Between1, [ 0.9290    0.6940    0.1250],'EdgeColor','none');
    % set(gca,'Layer','bottom');
    Between2 = [data(1:5,4)', fliplr(data(1:5,5)')];
    f2 = fill(x2, Between2, [0.3010    0.7450    0.9330],'EdgeColor','none');
    % set(gca,'Layer','bottom');
    % set(h,'edgecolor','white');
    h1 = get(gca,'Children');
    set(gca,'Children',[h1(3) h1(4) h1(1) h1(2)]);

    grid on;
    box on;
    xlim ([0.6 1.0]);
    ylim ([0 0.14]);



    title('\mu_s vs L_h , V = 0.2 [m/s]','fontsize',title_fontsize);
    xlabel('\mu_s' ,'fontsize',axis_fontsize);
    ylabel('L_h [m]','fontsize',axis_fontsize);

    set(gca,'Xtick',0.6:0.1:1.0,'Ytick',0:0.02:0.18,'fontsize',axis_fontsize);
    
    plot(1,0.1,'x','Color','r','linewidth',1.5);
    line([0 1],[0.1 0.1],'Color','r','LineStyle',':','linewidth',1.5);
    line([1 1],[0 0.1],'Color','r','LineStyle',':','linewidth',1.5);
    
    

    %%
    lgd = legend({'Wheel OP zone','leg OP zone','L_{h,max,wheeled}','L_{h,max,legged}','Stairs condition'},...
                'location','southeast');
    lgd.FontSize = legend_fontsize;
end
    %%


    % for visualization
    set(gcf,'Units','centimeters','position',[2 2 outputsize])%[5 5 2.5 1.8]

    fig.PaperPositionMode = 'auto';

    % for export
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 outputsize]);%[0 0 2.5 1.9],'PaperSize', [1.7 1.9]




%%

% %% plot sin critical points
% data = [
% 1	0.9	0.1	1.46 ;
% 1	0.9	0.2	0.77 ;
% 1	0.9	0.3	0.49 ;
% 1	0.9	0.4	0.38 ;
% 1	0.9	0.5	0.30 ;
% 1	0.9	0.6	0.25 ;
% 1	0.9	0.7	0.22 ;
% 1	0.9	0.8	0.20 ;
% 1	0.9	0.9	0.18 ;
% 1	0.9	1	0.18 ;
% 			
% 0.9	0.8	0.1	1.25 ;
% 0.9	0.8	0.2	0.68 ;
% 0.9	0.8	0.3	0.42 ;
% 0.9	0.8	0.4	0.32 ;
% 0.9	0.8	0.5	0.26 ;
% 0.9	0.8	0.6	0.22 ;
% 0.9	0.8	0.7	0.20 ;
% 0.9	0.8	0.8	0.18 ;
% 0.9	0.8	0.9	0.18 ;
% 0.9	0.8	1	0.16 ;
% 			
% 0.8	0.7	0.1	1.15 ;
% 0.8	0.7	0.2	0.60 ;
% 0.8	0.7	0.3	0.36 ;
% 0.8	0.7	0.4	0.28 ;
% 0.8	0.7	0.5	0.22 ;
% 0.8	0.7	0.6	0.20 ;
% 0.8	0.7	0.7	0.18 ;
% 0.8	0.7	0.8	0.16 ;
% 0.8	0.7	0.9	0.16 ;
% 0.8	0.7	1	0.14 ;
% 			
% 0.7	0.6	0.1	0.95 ;
% 0.7	0.6	0.2	0.48 ;
% 0.7	0.6	0.3	0.32 ;
% 0.7	0.6	0.4	0.24 ;
% 0.7	0.6	0.5	0.20 ;
% 0.7	0.6	0.6	0.18 ;
% 0.7	0.6	0.7	0.16 ;
% 0.7	0.6	0.8	0.14 ;
% 0.7	0.6	0.9	0.14 ;
% 0.7	0.6	1	0.12 ;
% 			
% 0.6	0.5	0.1	0.78 ;
% 0.6	0.5	0.2	0.40 ;
% 0.6	0.5	0.3	0.28 ;
% 0.6	0.5	0.4	0.20 ;
% 0.6	0.5	0.5	0.16 ;
% 0.6	0.5	0.6	0.14 ;
% 0.6	0.5	0.7	0.12 ;
% 0.6	0.5	0.8	0.12 ;
% 0.6	0.5	0.9	0.1 ;
% 0.6	0.5	1	0.1 ;
% 
% ];
%     
% 
% hold on;    
% for fric_index = 0.6:0.1:1
%     row_ind = data(:,1) == fric_index;  
%     plot(data(row_ind,3),data(row_ind,4),'*-');
% end
% 
% 
% 
% title('Critical points, Sin, V=0.4 [m/s]');
% xlabel('Frequency [Hz]');
% ylabel('Amplitude [m]');
% 
% xlim ([0 1.1]);
% % ylim ([0.05 0.20]);
% 
% 
% legend('Fs = 0.6, Fk = 0.5','Fs = 0.7, Fk = 0.6','Fs = 0.8, Fk = 0.7','Fs = 0.9, Fk = 0.8','Fs = 1.0, Fk = 0.9');

% %% different speed compare
% clear variables; clc; close all;
% 
% data = [
% 
% 0.4	1	0.9	0.1	1.46 ;
% 0.4	1	0.9	0.2	0.77 ;
% 0.4	1	0.9	0.3	0.49 ;
% 0.4	1	0.9	0.4	0.38 ;
% 0.4	1	0.9	0.5	0.30 ;
% 0.4	1	0.9	0.6	0.25 ;
% 0.4	1	0.9	0.7	0.22 ;
% 0.4	1	0.9	0.8	0.20 ;
% 0.4	1	0.9	0.9	0.18 ;
% 0.4	1	0.9	1	0.18 ;
% 				
% 0.2	1	0.9	0.1	1.44 ;
% 0.2	1	0.9	0.2	0.72 ;
% 0.2	1	0.9	0.3	0.52 ;
% 0.2	1	0.9	0.4	0.40 ;
% 0.2	1	0.9	0.5	0.32 ;
% 0.2	1	0.9	0.6	0.28 ;
% 0.2	1	0.9	0.7	0.24 ;
% 0.2	1	0.9	0.8	0.20 ;
% 0.2	1	0.9	0.9	0.18 ;
% 0.2	1	0.9	1	0.16 ;
% 				
% 				
% 0.4	0.6	0.5	0.1	0.78 ;
% 0.4	0.6	0.5	0.2	0.40 ;
% 0.4	0.6	0.5	0.3	0.28 ;
% 0.4	0.6	0.5	0.4	0.20 ;
% 0.4	0.6	0.5	0.5	0.16 ;
% 0.4	0.6	0.5	0.6	0.14 ;
% 0.4	0.6	0.5	0.7	0.12 ;
% 0.4	0.6	0.5	0.8	0.12 ;
% 0.4	0.6	0.5	0.9	0.1 ;
% 0.4	0.6	0.5	1	0.1 ;
% 				
% 0.2	0.6	0.5	0.1	0.80 ;
% 0.2	0.6	0.5	0.2	0.40 ;
% 0.2	0.6	0.5	0.3	0.28 ;
% 0.2	0.6	0.5	0.4	0.24 ;
% 0.2	0.6	0.5	0.5	0.20 ;
% 0.2	0.6	0.5	0.6	0.16 ;
% 0.2	0.6	0.5	0.7	0.14 ;
% 0.2	0.6	0.5	0.8	0.12 ;
% 0.2	0.6	0.5	0.9	0.12 ;
% 0.2	0.6	0.5	1	0.10 
% 
% ];
% 
% hold on;    
% 
% for vel_index = [0.2 0.4]
%     row_ind = (data(:,1) == vel_index) & (data(:,2) == 0.6);  
%     plot(data(row_ind,4),data(row_ind,5),'*-');
% end
% 
% 
% 
% title('Critical points, compare V=0.2 & 0.4 [m/s], Fs=0.6, Fk=0.5','FontSize',12);
% xlabel('Frequency [Hz]');
% ylabel('Amplitude [m]');
% 
% xlim ([0 1.1]);
% ylim ([0 1.6]);
% 
% 
% legend('Fs = 0.6, Fk = 0.5, V=0.2[m/s]','Fs = 0.6, Fk = 0.5, V=0.4[m/s]');
% %% test 
% figure('Color', 'w')
% box on
% x = 0.1:0.1:10;
% y = sin(x);
% area(x, y, 'FaceColor', 'b', 'EdgeColor', 'b')
% % set(gca,'Layer','top')


