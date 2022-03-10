%% Figure S1B and S1D
% - change in dendrite innervation vol. vs. L5PT dendrites
% - change in dendrite/axon innervation vol across cell types
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

%%
matlabPath = 'D:\udvary-et-al-2022\analysis\';
addpath(genpath([matlabPath 'functions\']));
tablePath = [matlabPath 'output\Tables\Robustness\']; 
figurePath = [matlabPath 'output\Figures\']; 

optStr = {'Dendrite','Axon'};
rShape = 0.2; 
rCV = 0.5; 
targetMean = 0; % No big difference between both options!
%targetMean = 1;

for opt = 1:2
    
    if opt==1
        % DendriteLength
        load([matlabPath 'data\robustness_morphologies\' ...
                            'DendriteAligned_Combination']);   
        numCellTypes = length(CellTypes);
    elseif opt==2
        % AxonLength
        load([matlabPath 'data\robustness_morphologies\Axon_Combination']);   
        numCellTypes = length(CellTypes);
    end
    
    shapeVal = [];
    ctshapeVal = [];
    CV_M_CT = nan(1,numCellTypes);
    CV_SD_CT = nan(1,numCellTypes);
    CV_M_CellTypes = nan(2,numCellTypes);
    CV_SD_CellTypes = nan(2,numCellTypes);

    sampleShapeThresholdList = nan(1,numCellTypes);
    M_ZeroVoxel = nan(1,numCellTypes);

    for ct = 1:numCellTypes
        
        if opt==1
            % [Voxels x Combinations]
            % -> Mean and SD over all samples 
            %   -> Mean: Samples 1 to n
            %   -> SD: Samples 2 to n
            sampleMean = LOOResults(ct).DendriteLengthMean; 
        elseif opt==2
            % [Voxels x Combinations]
            % -> Mean and SD over all samples 
            %   -> Mean: Samples 1 to n
            %   -> SD: Samples 2 to n
            sampleMean = LOOResults(ct).AxonLengthMean; 
        end
        idxValid = LOOResults(ct).idxValid;

        % Values with entire sample
        xTargetMean = sampleMean{end}; % Voxels x 1 (Combinations)
        idxInnervatedVoxels = xTargetMean>0;
        numVoxelsInnervated = sum(idxInnervatedVoxels);

        mTarget = zeros(size(idxValid));
        mTarget(idxValid) = xTargetMean; 
        mTarget = squeeze(sum(mTarget,2)); 

        x = find(sum(mTarget,1)>0);
        y = find(sum(mTarget,2)>0);

        xRange = [min(x)-2 max(x)+2];
        yRange = [min(y)-2 max(y)+2];

        mTarget = (mTarget==0);
        mTarget = mTarget(yRange(1):yRange(2),:);
        mTarget = mTarget(:,xRange(1):xRange(2));

        xZeroVoxels = []; 
        gZeroVoxels = [];
        sampleShapeThreshold = nan; 
        
        CV_M = nan(1,numSamples(ct)); 
        CV_SD = nan(1,numSamples(ct)); 

        M = nan(1,numSamples(ct));
        SD = nan(1,numSamples(ct)); 
        
        fileID = fopen([tablePath CellTypes{ct} '_' optStr{opt} '.csv'],'w');
        fprintf(fileID,'NumSamples,Mean,SD,CV,CV_Mean,CV_SD,');
        fprintf(fileID,'Min(InnervationVolume),Max(InnervationVolume),');
        fprintf(fileID,'Median(InnervationVolume),Mean(InnervationVolume),');
        fprintf(fileID,'SD(InnervationVolume),numInnervatedVoxels,\n');
        for i = 1:numSamples(ct)

            % Values for current sample
            x0 = sampleMean{i}; % Voxels x Combinations

            % Gives same results!
            if ~targetMean
                % Compute average always new for each sample
                xTargetMean = mean(x0,2); % Voxels x 1 (Combinations)
                idxInnervatedVoxels = xTargetMean>0;
                numVoxelsInnervated = sum(idxInnervatedVoxels);
            end

            % Extract relevant voxels
            x = x0(idxInnervatedVoxels,:);
            
            % Shape statistics
            xZeroVoxelsTMP = (numVoxelsInnervated-sum(x>0)) ...
                                    ./numVoxelsInnervated; 

            % Boxplot Shape [Combinations x 1]
            xZeroVoxels = [xZeroVoxels xZeroVoxelsTMP];
            gZeroVoxels = [gZeroVoxels i.*ones(size(xZeroVoxelsTMP))];
            
            % CV statistics           
            % Over all combinations
            m_tmp = mean(x,2)'; % Voxels x 1
            sd_tmp = std(x,[],2)'; % Voxels x 1
            
            CV_M(i) = mean(sd_tmp./m_tmp);
            CV_SD(i) = std(sd_tmp./m_tmp);
                               
            % Mean over Means
            M(i) = mean(x(:));
            % STD over Means
            sd_tmp2 = std(x,[],2); % Voxels x 1
            SD(i) = std(sd_tmp2,[],1);

            % Find threshold
            if isnan(sampleShapeThreshold) && median(xZeroVoxelsTMP)<rShape
                sampleShapeThreshold = i;      
            end
            
            if i==numSamples(ct)-1
                shapeVal = [shapeVal xZeroVoxelsTMP.*100];
                ctshapeVal = [ctshapeVal ct.*ones(size(xZeroVoxelsTMP))];
                
                CV_M_CT(ct) = CV_M(i);
                CV_SD_CT(ct) = CV_SD(i);
                
                M_ZeroVoxel(ct) = mean(xZeroVoxelsTMP.*100);
            end

            % Illustrate Example
            [~,idx] = min(abs(xZeroVoxelsTMP-median(xZeroVoxelsTMP)));
            x_example = x0(:,idx); 

            % Plot
            m = zeros(size(idxValid));
            m(idxValid) = x_example; 
            m = squeeze(sum(m,2)); 
            m = (m==0);
            m = m(yRange(1):yRange(2),:);
            m = m(:,xRange(1):xRange(2));
            
            % Print values
            valtmp = xZeroVoxelsTMP.*100; 
        
            % CV values
            fprintf(fileID,'%d,%.4f,%.4f,%.4f,',i,M(i),SD(i),SD(i)/M(i));
            fprintf(fileID,'%.4f,%.4f,',CV_M(i),CV_SD(i));
            % Innervation Values 
            fprintf(fileID,'%.4f,%.4f,',min(valtmp),max(valtmp));
            fprintf(fileID,'%.4f,%.4f,',median(valtmp),mean(valtmp));
            fprintf(fileID,'%.4f,%d,\n',std(valtmp),numVoxelsInnervated);
        end
        
        fclose('all');
        
        % Calculate ranges
        validRange = [sampleShapeThreshold:numSamples(ct)];
        if numel(validRange)==1
            validRange = [sampleShapeThreshold-1:numSamples(ct)];
        end
        
        CV_M_clean = CV_M(validRange(1:end-1));
        CV_SD_clean = CV_SD(validRange(1:end-1));
        idxCVThreshold = find(CV_M_clean<rCV,1);
        if isempty(idxCVThreshold)
            idxCVThreshold = nan; 
        end
        
        sampleCVThreshold = idxCVThreshold+validRange(1)-1;
        sampleShapeThresholdList(ct) = sampleShapeThreshold; 
        
        minInnervationVol = min(xZeroVoxels(gZeroVoxels==numSamples(ct)-1));
        minInnervationVolSample = min(xZeroVoxels(gZeroVoxels==validRange(1)));

        %% Figures
        f(ct) = figure(ct);   
        % Robustness of Shape
        subplot(2,2,opt*2-1);
        boxplot(xZeroVoxels.*100,gZeroVoxels,'Color','k','Symbol','k+');
        hold on;
        %plot([sampleShapeThreshold sampleShapeThreshold]-0.5,[0 100],'r--');
        set(gca,'Box','off','TickDir','out','XLim',[1 numSamples(ct)]-0.5, ...
            'YLim',[0 100]);
        xlabel(['#' CellTypes{ct} ' morphologies']);
        ylabel(['Change of ' optStr{opt} ' I-Volume [%]']);
        
        saveBoxPlotData([tablePath CellTypes{ct} '_' optStr{opt} ...
                    '_BoxPlotData.csv'],xZeroVoxels.*100,gZeroVoxels);

        % Robustness of relative distribution of shape -> CV of length
        subplot(2,2,opt*2);
        errorbar(validRange(1:end-1),CV_M(validRange(1:end-1)), ...
                    CV_SD(validRange(1:end-1)),'ko', ...
                    'MarkerFaceColor','k');
        set(gca,'Box','off','TickDir','out','XTick',validRange, ...
            'XLim',[validRange(1) validRange(end)]-0.5, ...
            'YLim',[0 max(CV_M(validRange)+CV_SD(validRange))*1.05]);
        xlabel(['#' CellTypes{ct} ' morphologies']);
        ylabel(['CV ' optStr{opt} ' density']);
        
        CVsufficient = 1; 
        if isnan(idxCVThreshold)
            idxCVThreshold = numel(CV_M_clean); 
            CVsufficient = 0; 
        end
        
        %% Print to terminal
        fprintf('>>> %s [%s] (n=%d)\n', ...
                    CellTypes{ct},optStr{opt},numSamples(ct));
        if sampleShapeThreshold==numSamples(ct)
            fprintf('    SAMPLE SIZE NOT SUFFICIENT FOR SHAPE GIVEN CRITERIA!\n');  
        end
        if CVsufficient==0
            fprintf('    SAMPLE SIZE NOT SUFFICIENT FOR CV GIVEN CRITERIA!\n');  
        end
        fprintf('    SamplesShape = %d (minInnervationVolume = %.2f) [minInnervationVolume = %.2f]\n', ...
                    sampleShapeThreshold,minInnervationVolSample,minInnervationVol);    
        fprintf('    SamplesVoxel = %d (CV = %.2f) [minCV = %.2f]\n', ...
                    sampleCVThreshold,CV_M_clean(idxCVThreshold),CV_M_clean(end)); 
                
        CV_M_CellTypes(1,ct) = CV_M_clean(idxCVThreshold);
        CV_M_CellTypes(2,ct) = CV_M_clean(end);
        CV_SD_CellTypes(1,ct) = CV_SD_clean(idxCVThreshold);
        CV_SD_CellTypes(2,ct) = CV_SD_clean(end);
    end
    
    % Result Figure Shape
    if opt==1
        ylimtmp = [0 20];
    else
        ylimtmp = [0 35];
    end
    
    if max(shapeVal(:))>ylimtmp(2)
        warning('%.2f larger than YLIM!',max(shapeVal(:))); 
    end
    
    fR = figure(12);
    subplot(2,2,opt*2-1);
    boxplot(shapeVal,ctshapeVal,'Color','k','Symbol','k+');
    hold on;
    plot([1 numCellTypes],[mean(M_ZeroVoxel) mean(M_ZeroVoxel)],'k-');
    set(gca,'Box','off','TickDir','out','XLim',[1 numCellTypes+1]-0.5, ...
        'YLim',ylimtmp,'XTick',[1:numCellTypes],'XTickLabel',CellTypes);
    xlabel(['Cell Types']);
    ylabel(['Change of ' optStr{opt} ' I-Volume [%]']);
    
    % Result Figure CV
    subplot(2,2,opt*2);
    errorbar([1:numCellTypes],CV_M_CT,CV_SD_CT,'ko', ...
                'MarkerFaceColor','k');
    hold on;
    plot([1 numCellTypes],[mean(CV_M_CT) mean(CV_M_CT)],'k-'); 
    set(gca,'Box','off','TickDir','out','XLim',[1 numCellTypes+1]-0.5, ...
        'YLim',[0 max(CV_M_CT+CV_SD_CT)*1.05],'XTick',[1:numCellTypes], ...
        'XTickLabel',CellTypes);
    xlabel(['Cell Types']);
    ylabel(['CV of ' optStr{opt} ' density']);
        
    % Save as table
    fileID = fopen([tablePath 'CellTypeSummary_' optStr{opt} '.csv'],'w');
    fprintf(fileID,'CellType,MaxNumSamples,Min(InnervationVolume),Max(InnervationVolume),');
    fprintf(fileID,'Median(InnervationVolume),Mean(InnervationVolume),');
    fprintf(fileID,'NumSamples(InnervationVolume<20%%),');
    fprintf(fileID,'CV_Mean(InnervationVolume<20%%),CV_SD(InnervationVolume<20%%1),');
    fprintf(fileID,'CV_Mean(maxNumSamples-1),CV_SD(maxNumSamples-1),\n');
    for ct = 1:numCellTypes
        valtmp = shapeVal(ctshapeVal==ct);
        valtmp1 = sampleShapeThresholdList(ct);
        valtmp2 = CV_M_CellTypes(1,ct);
        valtmp3 = CV_SD_CellTypes(1,ct); 
        if sampleShapeThresholdList(ct) == numSamples(ct)
            valtmp1 = nan; 
            valtmp2 = nan; 
            valtmp3 = nan; 
        end
        
        % Innervation volume
        fprintf(fileID,'%s,%d,%.4f,%.4f,',CellTypes{ct},numSamples(ct), ...
                                            min(valtmp),max(valtmp));
        fprintf(fileID,'%.4f,%.4f,',median(valtmp),mean(valtmp));
        fprintf(fileID,'%d,',valtmp1);
        
        % CV
        fprintf(fileID,'%.4f,%.4f,',valtmp2,valtmp3);
        fprintf(fileID,'%.4f,%.4f,\n',CV_M_CellTypes(2,ct),CV_SD_CellTypes(2,ct));
    end
    fclose('all');
    fprintf('%s: SamplingThreshold = %d\n',optStr{opt},numCombThreshold);
end

%%
for ct = 7 %L5PT  % 1:length(f)
    set(f(ct),'PaperPositionMode', 'auto','Position',[0 0 1200 300]);
    print(f(ct),'-dsvg','-r600',[figurePath CellTypes{ct} '_Summary.svg']); 
end

%%
set(fR,'PaperPositionMode', 'auto','Position',[0 0 800 500]);
print(fR,'-dsvg','-r600',[figurePath 'Morphologies_SummaryResult.svg']); 

function [] = saveBoxPlotData(filename,val,group)
    t = table(val',group','VariableNames',{'ChangeOfVolume','SampleSize'});
    writetable(t,filename);
end