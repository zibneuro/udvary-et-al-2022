%% Fig. 6C and Table S2
% Perin et al., 2011 (Fig 1E)
% L5tt -> L5tt Connection Probability drops with increasing Intersomatic 
% distance
% ISLimit_center = ([66.557 77.846 88.665 99.601 110.537 121.708 132.88 ...
%                             143.933 154.752]-61.009)./(155.692-61.009).*300;
% rangePerin = 25; 
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
load([matlabPath 'data\L5PT_L5PT_slice'], ...
            'CP_median','CP_25','CP_75','ISLimit','ISLimit_center');
        
%% Measurements
CP_measured = ([53.857 63.03 67.44 71.261 76.788 75.259 82.55 87.959 ...
                                81.257]-94.662)./(57.15-94.662).*0.2;
colExpmtl = [0 0.7 0];

%%
f1 = figure(1);
clf;
wdth = 0.3; 
lnwdth = 2;

fid = fopen([tablePath 'Perin.csv'],'w');
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
    'Ylim',[0 0.35],'TickDir','out','TickLength',[.02 .02],'Box','off', ...
    'XLim',[1-wdth*1.4 length(CP_median)+wdth*1.4], ...
    'XTick',[1:length(CP_median)],'XTickLabel',round(ISLimit_center));
set(f1,'PaperPositionMode', 'auto','Position',[0 0 500 400]);
print(f1,'-dsvg','-r600',[figurePath 'CP_L5PT_ISDist_Perin2011.svg']); 