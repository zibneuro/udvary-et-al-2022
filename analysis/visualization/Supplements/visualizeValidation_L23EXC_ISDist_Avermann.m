%% Script for Table S3
% Avermann et al.,2012
% L23->L23 Connection probability drops with increasing Intersomatic
% distance by -0.05 % per um
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
load([matlabPath 'data\L23EXC_L23EXC_slice.mat'], ...
            'CP_median','CP_25','CP_75','ISLimit','ISLimit_center');
        
%% Measurements
len = 54.008-[25.106 30.692 36.766 34.826];
CP_measured = len./54.008.*0.4;
colExpmtl = [0 0.7 0];

%%
f1 = figure(1);
clf;
wdth = 0.3; 
lnwdth = 2;

fid = fopen([tablePath 'Avermann.csv'],'w');
if fid==-1
    error('Cannot write file');
end
fprintf(fid,'IntersomaticDistance,p[measured],p[median],p[25th],p[75th],\n');

for i = 1:length(CP_median)
    
    plot([i-wdth i+wdth],[CP_median(i) CP_median(i)],'k-','LineWidth',lnwdth);
    hold on;
    rectangle('Position',[i-wdth CP_25(i) wdth*2 CP_75(i)-CP_25(i)],...
        'EdgeColor',[0 0 0],'LineWidth',lnwdth);

    % Measurement Plot
    plot(i,CP_measured(i),'o','MarkerFaceColor',colExpmtl, ...
        'MarkerSize',12,'MarkerEdgeColor','none');
    
    fprintf('Dist %.0f:\n Measured = %.2f\n Model: Med = %.2f 25thPrctl = %.2f 75thPrctl = %.2f\n', ...
            round(ISLimit_center(i)),CP_measured(i),CP_median(i), ...
            CP_25(i),CP_75(i));
        
    fprintf(fid,'%.0f,%.1f,%.1f,%.1f,%.1f,\n', ...
                round(ISLimit_center(i)),CP_measured(i)*100, ...
                CP_median(i)*100,CP_25(i)*100,CP_75(i)*100);
end
fclose(fid); 

set(gca,'YTick',[0:0.1:0.5], ...
    'Ylim',[0 0.5],'TickDir','out','TickLength',[.02 .02],'Box','off', ...
    'XLim',[1-wdth*1.4 length(CP_median)+wdth*1.4], ...
    'XTick',[1:length(CP_median)],'XTickLabel',round(ISLimit_center));
% set(f1,'PaperPositionMode', 'auto','Position',[0 0 200 300]);
% print(f1,'-dsvg','-r600',[figurePath 'CP_L23EXC_ISDist_Avermann2012.svg']); 