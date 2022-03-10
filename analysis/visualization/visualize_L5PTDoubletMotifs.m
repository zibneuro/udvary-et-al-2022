%% Figure 6E L5PT Doublet Motifs
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 

%%
load([matlabPath 'data\L5PT_Doublets.mat'], ...
        'doubleMotifDev','doubletMotifP','doubletMotifP_RND', ...
        'sliceIDs','doubletMotifLabel','numSamples');
    
fprintf('Samples = %d\n',sum(numSamples));

tblResult = table([1:numel(numSamples)]', ...
        doubletMotifP(:,1),doubletMotifP_RND(:,1),doubleMotifDev(:,1), ...
        doubletMotifP(:,2),doubletMotifP_RND(:,2),doubleMotifDev(:,2), ...
        doubletMotifP(:,3),doubletMotifP_RND(:,3),doubleMotifDev(:,3), ...
        numSamples','VariableNames',{'SlideID', ...
        'p(model_unidirection)','p(random_unidirection)','dev(unidirection)', ...
        'p(model_bidirection)','p(random_bidirection)','dev(bidirection)', ...
        'p(model_unconnected)','p(random_unconnected)','dev(unconnected)', ...
        'Samples'});
writetable(tblResult,[tablePath 'L5PT_Doublets.csv']);
    
% Compare to Song et al., 2005
pavg_Song = 931/8050;
% Frequency of unidirectional, bidirectional, unconnected motifs
freqMotifSong = [495 218 3312]; 
pMotifSong = [2*pavg_Song*(1-pavg_Song) pavg_Song^2 (1-pavg_Song)^2 ];
devSong = freqMotifSong./sum(freqMotifSong)./pMotifSong;
devSong_err = [1/1.116523 1.076336 1.05].*devSong; % measured from paper
devModel_avg = mean(doubleMotifDev,1);
devModel_min = min(doubleMotifDev);
devModel_max = max(doubleMotifDev);

% Log10 and normalized
devSongLog = log10(devSong);
devSongLog_err = log10(devSong_err); 
devModelLog = mean(log10(doubleMotifDev),1); 
devSongLog_err = devSongLog_err./max(devSongLog); 
devSongLog = devSongLog./max(devSongLog); 
devModelLog_min = log10(devModel_min)./max(devModelLog);
devModelLog_max = log10(devModel_max)./max(devModelLog);
devModelLog = devModelLog./max(devModelLog); 

tblMeasuredDoublet = table();
tblMeasuredDoublet.Motif = {'1-Connection','2-Connections','0-Connections'}';
tblMeasuredDoublet.MotifCount = freqMotifSong';
tblMeasuredDoublet.MotifProbability = freqMotifSong'./sum(freqMotifSong);
tblMeasuredDoublet.MotifProbabilityRandom = pMotifSong';
tblMeasuredDoublet.Deviation = devSong'; 
writetable(tblMeasuredDoublet,[tablePath 'SongEtAl_DoubletMeasurement.csv']);

%% Figure 
f1 = figure(1);
clf;
barLen = 0.1;
colGreen = [20 205 20]./255; 

% Raw data
subplot(2,1,1);
h = bar([devSong; devModel_avg]',1,'Basevalue',1,'EdgeColor','none');
h(1).FaceColor = colGreen;
h(2).FaceColor = 'k';
hold on;
% plot errorbars
% Model: min/max over 20 slices
% measured: bootstrap errorbar
for i = 1:3
    % Model
    plot([i i]+0.15,[devModel_min(i) devModel_max(i)],'-', ...
            'Color',[0.5 0.5 0.5]);
    plot([-1 1].*barLen+i+0.15,[devModel_max(i) devModel_max(i)], ...
            '-','Color',[0.5 0.5 0.5]);
    plot([-1 1].*barLen+i+0.15,[devModel_min(i) devModel_min(i)], ...
            '-','Color',[0.5 0.5 0.5]);
    % Measured
    plot([i i]-0.15,[devSong(i) devSong_err(i)],'-', ...
            'Color',colGreen);
    plot([-1 1].*barLen+i-0.15,[devSong_err(i) devSong_err(i)], ...
            '-','Color',colGreen);
end

set(gca,'Box','off','TickDir','out', ...
    'XTick',[1:3],'XTickLabel',doubletMotifLabel, ...
    'YTick',[0:1:4],'YScale','log');
ylabel('deviation');

% Log10 normalized to max
subplot(2,1,2);
h = bar([devSongLog; devModelLog]',1,'Basevalue',0,'EdgeColor','none');
h(1).FaceColor = colGreen;
h(2).FaceColor = 'k';
hold on;
for i = 1:3
    % Model
    plot([i i]+0.15,[devModelLog_min(i) devModelLog_max(i)],'-', ...
            'Color',[0.5 0.5 0.5]);
    plot([-1 1].*barLen+i+0.15,[devModelLog_max(i) devModelLog_max(i)], ...
            '-','Color',[0.5 0.5 0.5]);
    plot([-1 1].*barLen+i+0.15,[devModelLog_min(i) devModelLog_min(i)], ...
            '-','Color',[0.5 0.5 0.5]);
    % Measured
    plot([i i]-0.15,[devSongLog(i) devSongLog_err(i)],'-', ...
            'Color',colGreen);
    plot([-1 1].*barLen+i-0.15,[devSongLog_err(i) devSongLog_err(i)], ...
            '-','Color',colGreen);
end
set(gca,'Box','off','TickDir','out', ...
    'XTick',[1:3],'XTickLabel',doubletMotifLabel, ...
    'YTick',[-1:0.5:1],'YScale','linear');
ylabel('deviation');

set(f1,'PaperPositionMode', 'auto','Position',[0 0 200 400]);
print(f1,'-dsvg','-r600',[figurePath 'L5tt_TwoNeuronMotifSpectra.svg']);