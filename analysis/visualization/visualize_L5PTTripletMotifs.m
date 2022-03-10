%% Figure 6E
% Visualize L5PT triplet motif spectra
% normalized by doublet motif spectra
% over 20 slices
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
newIdxOrder = [1 3 5 6 4 7 2 8 10 11 12 13 9 14 15 16]; 

%% Song et al
% Measured counts of motifs
idxMeasuredSort = [16 15 14 12 13 10 11 9 8 7 5 4 6 3 2 1]; 
measuredMotifCount = [1375 579 274 33 25 41 41 24 17 9 4 6 4 5 6 3];
measuredMotifCount = measuredMotifCount(idxMeasuredSort); 
measuredMotifP = measuredMotifCount./sum(measuredMotifCount); 

% 1. Approach [Normalize by doublet motifs]
% Use Song et al approach for their measured data and our approach
% for model's prediction (delta Probability)
% Song et al.: random p is normalized by probabiltiy of
% bidirectional and unidirectional edge probability
% where bidirectional is not unidirectional^2 but what they measured!
% Frequency of unidirectional, bidirectional, unconnected doublet motifs
% from page 0510:
freqMotifSong = [495 218 3312]; 
m = freqMotifSong./sum(freqMotifSong); % respective probability per motif
m(1) = m(1)./2; % normalize unidirectional motif (both ways)
% For reference: probability numbers from paper (match manual calculation)
% m(1) = 0.0615; % unidirectinal (from paper, Fig 4A)
% m(2) = 0.0542; % bidirectional (from paper, Fig 4A)
% m(3) = 1-m(1)*2-m(2); % no connection 
measuredPRandom = computeMotifProbabilitiesBasedOnDoubletMotifs(m);
devSong = measuredMotifP./measuredPRandom;           
devSong = devSong';
scaleFac = [1.653298391 1.45 1.45 1.354960553 ...
        1.5 1.4 1.45 1.15 1/1.1 ...
        1.05 1.05 1.1 1/1.3 1.05 1/1.05 1]'; 
devSong_err = scaleFac.*devSong; % measured from paper

devSongLog = log10(devSong);
devSongLog_err = log10(devSong_err); 
devSongLog_err = devSongLog_err./max(devSongLog); 
devSongLog = devSongLog./max(devSongLog); 

% Save measured data in table
tblMeasuredTriplet = table();
tblMeasuredTriplet.MotifID = [1:numel(measuredMotifCount)]'; 
tblMeasuredTriplet.MotifCount = measuredMotifCount';
tblMeasuredTriplet.MotifProbability = measuredMotifP';
tblMeasuredTriplet.MotifProbabilityFromDoublets = measuredPRandom';
tblMeasuredTriplet.Deviation = devSong;
writetable(tblMeasuredTriplet,[tablePath 'SongEtAl_TripletMeasurement.csv']);

%% Model
load([matlabPath 'data\L5PT_Doublets.mat'],'doubletMotifP');
load([matlabPath 'data\L5PT_Triplets.mat'],'p_avg','pMotif64');
numSlices = length(p_avg); 
pModel = collapse64MotifTo16Motif(pMotif64');
pModel = pModel(newIdxOrder,:); 

% 1. Approach: Normalize by DoubletMotifSpectra
devModel = nan(16,numSlices);
pModel_Rnd = nan(16,numSlices); 
for i = 1:numSlices   
    pRnd_tmp = computeMotifProbabilitiesBasedOnDoubletMotifs( ...
                    doubletMotifP(i,:)./[2 1 1]);
    devModel(:,i) = pModel(:,i)./pRnd_tmp';    
    pModel_Rnd(:,i) = pRnd_tmp';   
end

devModel_avg = mean(devModel,2);
devModel_min = min(devModel,[],2);
devModel_max = max(devModel,[],2);

devModelLog = mean(log10(devModel),2); 
devModelLog_min = log10(devModel_min)./max(devModelLog);
devModelLog_max = log10(devModel_max)./max(devModelLog);
devModelLog = devModelLog./max(devModelLog); 

tblResult = table([1:size(devModel,2)]',pModel(1:15,:)',pModel_Rnd(1:15,:)', ...
        devModel(1:15,:)', ...
        'VariableNames',{'SlideID','p(model)_motif', ...
        'p(random_doublet)_motif','deviation_motif'});
writetable(tblResult,[tablePath 'L5PT_Triplets.csv']);

%% Figure
f1 = figure(1);
clf;
barLen = 0.1;
colGreen = [20 205 20]./255; 

subplot(211);
b = bar([1:16],[devSong devModel_avg],...
    1,'EdgeColor','none','FaceColor','k','BaseValue',1); 
b(1).FaceColor = colGreen;
b(2).FaceColor = 'k';
hold on;
% plot errorbars
% Model: min/max over 20 slices
% measured: bootstrap errorbar
for i = 1:16
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
set(gca,'TickDir','out','Box','off','XLim',[0.5 15.5],'XTick', ...
            [1:16],'YScale','linear'); 
ylabel('deviation');

subplot(212);
b = bar([1:16],[devSongLog devModelLog],...
    1,'EdgeColor','none','FaceColor','k','BaseValue',0); 
b(1).FaceColor = colGreen;
b(2).FaceColor = 'k';
hold on;
for i = 1:16
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
set(gca,'TickDir','out','Box','off','XLim',[0.5 15.5],'XTick', ...
            [1:16],'YScale','linear'); 
ylabel('deviation (a.u.)');

set(f1,'PaperPositionMode','auto','Position',[0 0 500 350]);
print(f1,'-dsvg','-r600',[figurePath 'L5ttMotifs.svg']);