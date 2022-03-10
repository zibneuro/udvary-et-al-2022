%% Figure 6D
% measured vs predicted connection probability across 89 measurements
% as scatter/correlation plot
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\'];  

%%
tbl = readtable([matlabPath 'data\cellular_connectivity\stats\summary.csv']); 
CP_measured = tbl.empirical';
CP_AVG_predicted = tbl.avg_slices_merged';
CP_SD_predicted = tbl.std_slices_merged';
deviationValues = (CP_measured-CP_AVG_predicted)./CP_SD_predicted;

[r,p] = corr(CP_AVG_predicted',CP_measured');
b = CP_AVG_predicted'\CP_measured'; % Simple Linear regression

% Compute 95% Prediction Interval for new observations
% (might not be meaningful in this context here)
mdl = fitlm(CP_AVG_predicted',CP_measured','linear','intercept',false);
[~,yci] = predict(mdl,CP_AVG_predicted','Prediction','observation');
fprintf('r = %.4f p = %.4e (n = %d)\n',r,p,numel(CP_AVG_predicted));
fprintf('Devation Values: M = %.2f SD = %.2f Range = [%.2f %.2f]\n', ...
         mean(abs(deviationValues)),std(abs(deviationValues)), ...
         min(abs(deviationValues)),max(abs(deviationValues)));
fprintf('Deviation<0.5: %d (%.2f)\n',sum(abs(deviationValues)<0.5), ...
    sum(abs(deviationValues)<0.5)/numel(deviationValues));
fprintf('Deviation<1: %d (%.2f)\n',sum(abs(deviationValues)<1), ...
    sum(abs(deviationValues)<1)/numel(deviationValues));   

fprintf('#measurements = %d (#studies = %d)\n',numel(CP_AVG_predicted), ...
            numel(unique(tbl.author_year)));
     
%% Figures Correlation
f1 = figure(1);
clf;
% maxVal = round(max([CP_AVG_predicted CP_measured]+0.05),1);
maxVal = 0.7; 

if maxVal<max([CP_AVG_predicted CP_measured])
   warning('Axis limit cuts out values'); 
end

plot(CP_AVG_predicted,CP_measured,'o','MarkerFaceColor','k', ...
            'MarkerEdgeColor','none');
hold on;
plot([0 maxVal],[0 maxVal],'k--');
plot(CP_AVG_predicted,CP_AVG_predicted.*b,'r-');
plot(CP_AVG_predicted,yci,'r-');

% Examples for paper
% 3: VPM -> L4 Bruno & Sakman 2006 (red)
% 46: L4sp->L4sp Sun et al. ,2006 (anatomy EM) (green)
% 36: L4->L23 (Yoshimura, Optical Stimulation)  (blue)
% 65: L5tt-L5tt (Perin, Berger et al. 2011,)    (cyan)
IDList = [3 46 36 65];
col = 'rgbc'; 

for i = 1:length(IDList)
    plot(CP_AVG_predicted(IDList(i)),CP_measured(IDList(i)),'o', ...
        'MarkerFaceColor','none','MarkerEdgeColor',col(i));
end

set(gca,'YLim',[-0.03 maxVal],'XLim',[-0.03 maxVal],'Box','off', ...
    'TickDir','out','XTick',[0:0.1:1],'YTick',[0:0.1:1]);
xlabel('predicted p');
ylabel('measured p');
axis square;

set(f1,'PaperPositionMode', 'auto','Position',[0 0 300 300]);
print(f1,'-dsvg','-r600',[figurePath 'MeasuredVsPredictedCorrPlot.svg']);