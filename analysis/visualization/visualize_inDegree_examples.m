%% Figure 4E
% Correlation (L5ACC,L5PT)->L5PT
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figurePath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\']; 
addpath(genpath([matlabPath 'functions\']));
load([matlabPath 'data\CellTable_RBC_20.mat'],'CellTable');
idxDel = strcmp(CellTable.nearest_column,'C2');
CellTable(~idxDel,:)= []; 
idxSep = strcmp(CellTable.region,'S1_Septum_C2');
idxL5 = (CellTable.layer==5);

%% Example Figure
cellType1List = {'L5PT','L5IT','L5IT',};    % (A) [Presynapse]
cellType2List = {'L6ACC','L6ACC','L6ACC',};    % (B) [Presynapse]
cellType3List = {'L5PT','L6ACC','L5IT'};    % (C) [Postsynapse]

f1 = figure(1);
clf;

for i = 1
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
    idxtmp = ismember(I1.I.PostCellID,CellTable.neuronID(idxSep));      

    % C2
    syn1 = sum(I1.I.I,1);
    syn2 = sum(I2.I.I,1);

    % C2 Septum
    syn1Sep = syn1(idxtmp);
    syn2Sep = syn2(idxtmp);
    
    % C2 inside column
    syn1BC = syn1(~idxtmp);
    syn2BC = syn2(~idxtmp);    
    
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
    fprintf(fid,'\ninDegree(%s_%s[Septum]),inDegree(%s_%s[Septum]),\n', ...
                    preType1,postType,preType2,postType);
    for ii = 1:numel(syn1Sep)
        fprintf(fid,'%.4f,%.4f,\n',syn1Sep(ii),syn2Sep(ii));
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
    
    % Corr Septal
    [r,p,RL,RU] = corrcoef(syn1Sep,syn2Sep,'Alpha',0.05);
    if p(1,2)<0.05
        str = sprintf('r = %.2f [%.2f %.2f] (*) p = %.2e (n=%d)', ...
            r(1,2),RL(1,2),RU(1,2),p(1,2),numel(syn1Sep)); 
    else
        str = sprintf('r = %.2f [%.2f %.2f] p = %.2e (n=%d)', ...
            r(1,2),RL(1,2),RU(1,2),p(1,2),numel(syn1Sep)); 
    end
    fprintf('(%s,%s)->%s[Septum]: %s\n',preType1,preType2,postType,str); 
    
    % Corr Inside BC
    [r,p,RL,RU] = corrcoef(syn1BC,syn2BC,'Alpha',0.05);
    if p(1,2)<0.05
        str = sprintf('r = %.2f [%.2f %.2f] (*) p = %.2e (n=%d)', ...
            r(1,2),RL(1,2),RU(1,2),p(1,2),numel(syn1BC)); 
    else
        str = sprintf('r = %.2f [%.2f %.2f] p = %.2e (n=%d)', ...
            r(1,2),RL(1,2),RU(1,2),p(1,2),numel(syn1BC)); 
    end
    fprintf('(%s,%s)->%s[BC]: %s\n',preType1,preType2,postType,str); 

    % Polyfit
    p = polyfit(syn1,syn2,1);
    yfit = polyval(p,[min(syn1) max(syn1)]); 
    % Polyfit
    p = polyfit(syn1Sep,syn2Sep,1);
    yfitSep = polyval(p,[min(syn1Sep) max(syn1Sep)]); 
    % Polyfit
    p = polyfit(syn1BC,syn2BC,1);
    yfitBC = polyval(p,[min(syn1BC) max(syn1BC)]); 
    
    X = [syn1 syn2];
    mx = max(X(:)); 
    mx_x = max(syn1);
    mx_y = max(syn2); 
    
    % Plot
    subplot(3,1,i);     
    if i~=3 
        scatter(syn1,syn2,10,'filled','MarkerEdgeColor','none', ...
                        'MarkerFaceColor',[0.35 0.35 0.35],'MarkerFaceAlpha',0.3);
        hold on;
        plot([min(syn1) max(syn1)],yfit,'-','LineWidth',1,'Color',col); 
   else % In Layer 5 if i==3
        scatter(syn1BC,syn2BC,10,'filled','MarkerEdgeColor','none', ...
                        'MarkerFaceColor',col,'MarkerFaceAlpha',0.3);
        hold on;
        scatter(syn1Sep,syn2Sep,10,'filled','MarkerEdgeColor','none', ...
                        'MarkerFaceColor',[139 197 63]./255,'MarkerFaceAlpha',0.3);
        plot([min(syn1BC) max(syn1BC)],yfitBC,'-','LineWidth',1, ...
                        'Color',col); 
        plot([min(syn1Sep) max(syn1Sep)],yfitSep,'-','LineWidth',1, ...
                        'Color',[139 197 63]./255); 
   end
    XRange = [-0.03*mx_x mx_x*1.1]; 
    YRange = [-0.03*mx_y mx_y*1.1]; 
    plot([0 mx],[0 mx],':','Color',[0.5 0.5 0.5]);
    
    xlabel(['in-degree ' preType1 '[C2]']); 
    ylabel(['in-degree ' preType2 '[C2]']);
    title(['onto ' postType '[C2]']);
    set(gca,'TickDir','out','Box','off','XLim',XRange, ...
        'YLim',YRange); %'XTick',Ticks(i,:),'YTick',Ticks(i,:));
end

set(f1,'PaperPositionMode', 'auto','Position',[0 0 300 700]);
print(f1,'-dsvg','-r600',[figurePath 'Correlation_Examples_Subsamples.svg']);

fclose('all'); 