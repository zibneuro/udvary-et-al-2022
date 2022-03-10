%% Extract Table of predicted vs empirical connection probabilities
% Table S1
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
tablePath = [matlabPath 'output\Tables\']; 
dataPath = [matlabPath 'data\cellular_connectivity\stats\'];

% From old manuscript
sortIDX = [3 1 4 2 5 6 7 8 9 10 11 12 13 14 15 16 22 20 18 23 17 19 21 ...
    24 25 26 27 28 29 30 31 32 33 34 35 36 38 37 39 40 41 42 43 44 45 ...
    46 47 48 49 50 54 51 53 52 55 56 57 58 59 60 61 62 64 66 65 63 67 ...
    68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 87 86 88 89];
numSamples = [40 62 14 21 24 11 18 9 11 950 878 183 208 211 108 50 95 ...
    nan 542 112 247 112 235 774 98 51 182 513 170 29 87 164 61 208 64 ...
    50 172 25 89 1046 nan 276 136 93 528 24 94 146 24 29 150 500 1450 ...
    163 98 118 66 86 36 51 86 36 225 12 3446 8050 209 89 275 934 175 ...
    158 104 167 137 174 269 555 100 50 64 94 160 100 102 532 27 43 40]; 
            
% Final sorting
species = {'rat' 'rat' 'rat' 'rat' 'rat' 'rat' 'rat' 'rat' 'rat' 'mouse' ...
    'mouse' 'mouse' 'mouse' 'mouse' 'mouse' 'mouse' 'mouse' 'rat' 'rat' ...
    'mouse' 'rat' 'rat' 'mouse' 'mouse' 'mouse' 'mouse' 'mouse' 'mouse' ...
    'mouse' 'rat' 'mouse' 'mouse' 'mouse' 'mouse' 'rat' 'rat' 'mouse' ...
    'rat' 'rat' 'mouse' 'rat' 'mouse' 'mouse' 'mouse' 'rat' 'rat' 'rat' ...
    'rat' 'rat' 'rat' 'mouse' 'rat' 'rat' 'rat' 'mouse' 'mouse' 'mouse' ...
    'mouse' 'mouse' 'mouse' 'mouse' 'mouse' 'mouse' 'mouse' 'rat' 'rat' ...
    'mouse' 'mouse' 'mouse' 'mouse' 'mouse' 'mouse' 'mouse' 'mouse' ...
    'mouse' 'mouse' 'mouse' 'mouse' 'mouse' 'mouse' 'mouse' 'mouse' ...
    'mouse' 'mouse' 'rat' 'mouse' 'rat' 'mouse' 'mouse'};
age = {'24-35' 'adult' '24-35' 'adult' '24-35' '24-35' 'adult' 'adult' ...
    'adult' '18-21' '18-30' '18-21' '18-21' '18-21' '18-21' '18-21' ...
    '17-22' '17-23' '14-16' '60+' 'adult' '21-26' 'adult' '21-30' '60+' ...
    '60+' '18-21' '18-21' '18-21' 'adult' '18-21' '18-21' '18-21' '18-21' ...
    '17-23' '21-26' '18-21' 'adult' '14-21' '18-21' '17-23' '18-21' ...
    '18-21' '18-21' 'adult' '21-35' '12-15' '13-15' '21-35' 'adult' ...
    '60+' '14-16' '14-16' 'adult' '60+' '14-17' '60+' '14-17' '60+' ...
    '60+' '14-17' '60+' '14-17' '60+' '14-16' '12-20' '18-21' '18-21' ...
    '18-21' '18-21' '18-21' '18-21' '18-21' '18-21' '18-21' '18-21' ...
    '14-37' '18-21' '18-21' '18-21' '18-21' '18-21' '18-21' '18-21' ...
    '14-21' '18-21' 'adult' '20-33' '20-33'};

%%
tbl = readtable([dataPath 'summary.csv']); 
CP_measured = tbl.empirical';
CP_AVG_predicted = tbl.avg_slices_merged';
CP_SD_predicted = tbl.std_slices_merged';
deviationValues = (CP_measured-CP_AVG_predicted)./CP_SD_predicted;
N = size(tbl,1);
pBins = [0:1:100]'; % 1% bins

% Print table with following columns
% Create new table
tblFinal = table();
tblFinal.references = tbl.author_year;
tblFinal.preLayer = repmat({''},N,1); 
tblFinal.preType = repmat({''},N,1); 
tblFinal.postLayer = repmat({''},N,1); 
tblFinal.postType = repmat({''},N,1); 
tblFinal.species = repmat({''},N,1); 
tblFinal.age = repmat({''},N,1); 
tblFinal.invivo = repmat({'no'},N,1); 
tblFinal.numSamples = nan(N,1); 

% Connection probability values
tblFinal.p_measured = round(tbl.empirical*100);
tblFinal.p_avg = round(tbl.avg_slices_merged*100);
tblFinal.p_sd = round(tbl.std_slices_merged*100);

% Validation
tblFinal.dev = round(deviationValues'.*100)./100; 
tblFinal.prctl = nan(N,1); 

m = load([matlabPath 'data\MeasuredCP.mat']);
x.numSamples = m.NumSamples;
x.cp = ((m.CP_max+m.CP_min)./2).*100; 
x.year = m.Year;
x.authors = m.Authors; 
x.species = m.Species;
x.age = m.AnimalAge; 
clear m; 

%%
for i = 1:N
    if contains(tbl.descriptor(i),'invivo')
       tblFinal.invivo(i) = {'yes'}; 
    end
    splitStr = split(tbl.descriptor{i},["-","_"]);
    tblFinal.preLayer(i) = splitStr(2);
    tblFinal.preType(i) = splitStr(3);
    tblFinal.postLayer(i) = splitStr(4);
    tblFinal.postType(i) = splitStr(5);
    if contains(tbl.descriptor(i),'septum')
       tblFinal.postLayer(i) = {[splitStr{4} 'Septum']}; 
    end
    
    % Get percentile
    filename = [dataPath tbl.descriptor{i} '_histogram.txt'];
    tbltmp = readtable(filename); 
    n1 = sum(tbltmp.Var1(pBins<=tblFinal.p_measured(i)));
    n2 = sum(tbltmp.Var1);
    tblFinal.prctl(i) = round(n1/n2*100); 
end

tblFinal.preType = strrep(tblFinal.preType,'EXC','');
tblFinal.postType = strrep(tblFinal.postType,'EXC','');
tblFinal.preType = strrep(tblFinal.preType,'VPM','');

tblFinal = tblFinal(sortIDX,:);
tblFinal.numSamples = numSamples'; 
tblFinal.age = age';
tblFinal.species = species';

writetable(tblFinal,[tablePath 'empirical_vs_predicted_cp.csv']); 

%%
[r,p] =corr(tblFinal.p_avg,tblFinal.p_measured);
fprintf('All: R = %.2f p = %.2e (n=%d)\n',r,p,numel(tblFinal.p_avg));

idxRat = strcmp(tblFinal.species,'rat');
idxMouse = strcmp(tblFinal.species,'mouse'); 

[r,p] =corr(tblFinal.p_avg(idxRat),tblFinal.p_measured(idxRat));
fprintf('Rat Only: R = %.2f p = %.2e (n=%d)\n',r,p,sum(idxRat));

[r,p] =corr(tblFinal.p_avg(idxMouse),tblFinal.p_measured(idxMouse));
fprintf('Mouse Only: R = %.2f p = %.2e (n=%d)\n',r,p,sum(idxMouse));