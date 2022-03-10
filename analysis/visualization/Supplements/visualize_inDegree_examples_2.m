%% Figure S6C
% - InDegrees (L6ACC,L5IT)->L5IT
% - InDegrees (L6ACC,L5IT)->L4sp
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
addpath(genpath([matlabPath 'functions\']));

%% Example Figure
cellType1List = {'L5IT','L5IT'};    % (A) [Presynapse]
cellType2List = {'L6ACC','L6ACC'};    % (B) [Presynapse]
cellType3List = {'L5IT','L4sp'};    % (C) [Postsynapse]

f1 = figure(1);
clf;

for i = 1:length(cellType1List)

    preType1 = cellType1List{i}; % A (Presynapse)
    preType2 = cellType2List{i}; % B (Presynapse)
    postType = cellType3List{i}; % C (Postsynapse)
    col = [0 0 0];
    
    I1 = load([matlabPath 'data\CellMatrix_C2\' ...
                preType1 '_' postType '.mat']);
    I2 = load([matlabPath 'data\CellMatrix_C2\' ...
                preType2 '_' postType '.mat']);
    if sum(I1.I.PostCellID==I2.I.PostCellID)~=numel(I2.I.PostCellID)
        error('PostCellIDs do not match!');
    end    

    % C2
    syn1 = sum(I1.I.I,1);
    syn2 = sum(I2.I.I,1);

    % Save raw data
    fid = fopen([tablePath 'InDegree_' preType1 '_' preType2 '_to_' ...
                    postType '.csv'], 'w+');
    if fid==-1
       error('Cannot write .csv file'); 
    end
    fprintf(fid,'inDegree(%s_%s),inDegree(%s_%s),\n', ...
                    preType1,postType,preType2,postType);
    for ii = 1:numel(syn1)
        fprintf(fid,'%.4f,%.4f,\n',syn1(ii),syn2(ii));
    end
    fclose(fid);    
    
    % Corr
    [r,p,RL,RU] = corrcoef(syn1,syn2,'Alpha',0.05);
    
    if p(1,2)<0.05
        str = sprintf('r = %.2f [%.2f %.2f] (*) p = %.2e (n=%d)', ...
            r(1,2),RL(1,2),RU(1,2),p(1,2),numel(syn1)); 
    else
        str = sprintf('r = %.2f [%.2f %.2f] p = %.2e (n=%d)', ...
            r(1,2),RL(1,2),RU(1,2),p(1,2),numel(syn1)); 
    end
    fprintf('(%s,%s)->%s: %s\n',preType1,preType2,postType,str); 
    
    % Polyfit
    p = polyfit(syn1,syn2,1);
    yfit = polyval(p,[min(syn1) max(syn1)]); 

    X = [syn1 syn2];
    mx = max(X(:)); 
    mx_x = max(syn1);
    mx_y = max(syn2); 
    
    % Plot
    subplot(2,1,i);     
    scatter(syn1,syn2,10,'filled','MarkerEdgeColor','none', ...
                    'MarkerFaceColor',[0.35 0.35 0.35],'MarkerFaceAlpha',0.3);
    hold on;
    plot([min(syn1) max(syn1)],yfit,'-','LineWidth',1,'Color',col); 

    XRange = [-0.03*mx_x mx_x*1.1]; 
    YRange = [-0.03*mx_y mx_y*1.1]; 
    plot([0 mx],[0 mx],':','Color',[0.5 0.5 0.5]);
    
    xlabel(['in-degree ' preType1 '[C2]']); 
    ylabel(['in-degree ' preType2 '[C2]']);
    title(['onto ' postType '[C2]']);
    set(gca,'TickDir','out','Box','off','XLim',XRange, ...
        'YLim',YRange); %'XTick',Ticks(i,:),'YTick',Ticks(i,:));
end

set(f1,'PaperPositionMode', 'auto','Position',[0 0 300 400]);
print(f1,'-dsvg','-r600',[figurePath 'Correlation_Examples_Supplements.svg']);