%% Figure S4
% Example of Motif Deviations for different groups of neurons
% - different intersomatic distances
% - different neuron populations in L4 of C2
% - different neuron populations in L2/3 of C2
% - different neuron populations across different columns
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
load([matlabPath 'data\TripletMotifs_All.mat'],'tbl'); 
figurePath = [matlabPath 'output\Figures\']; 

%%
f1 = figure(1);
clf; 

% Intersomatic distances within C2
% 0 to 100 (389)
% 100 to 200 (390)
% 200 to inf (391)
subplot(2,2,1);
idx = [389 390 391];
dev = tbl.pMotif(idx,:)./tbl.pMotifRnd(idx,:);
for i = 1:3
    fprintf('%s %s %s\n',tbl.A{idx(i)},tbl.B{idx(i)},tbl.C{idx(i)});
end
fprintf('-----\n');
h = bar(1:15,[dev(:,1:15)]',1,'BaseValue',1,'EdgeColor','none'); 
h(1).FaceColor = 0.1.*ones(1,3);
h(2).FaceColor = 0.3.*ones(1,3);
h(3).FaceColor = 0.7.*ones(1,3);
set(gca,'TickDir','out','Box','off','XLim',[0.5 15.5],'XTick', ...
        [1:16],'YScale','log','Ylim',[min(dev(:)) max(dev(:))]); 
title('Intersomatic distances');

% Layer 4
% L4EXC (1804)
% L4INH (849)
% L4SP (44)
subplot(2,2,4);
idx = [1804 849 44];
dev = tbl.pMotif(idx,:)./tbl.pMotifRnd(idx,:);
for i = 1:3
    fprintf('%s %s %s\n',tbl.A{idx(i)},tbl.B{idx(i)},tbl.C{idx(i)});
end
fprintf('-----\n');
h = bar(1:15,[dev(:,1:15)]',1,'BaseValue',1,'EdgeColor','none'); 
h(1).FaceColor = 0.1.*ones(1,3);
h(2).FaceColor = 0.3.*ones(1,3);
h(3).FaceColor = 0.7.*ones(1,3);
set(gca,'TickDir','out','Box','off','XLim',[0.5 15.5],'XTick', ...
        [1:16],'YScale','log','Ylim',[min(dev(:)) max(dev(:))]);  
title('Layer 4');

% Layer 2/3
% L2EXC (1363)
% L3EXC (1573)
% L3INH (713)
subplot(2,2,2);
idx = [1363 1573 713];
dev = tbl.pMotif(idx,:)./tbl.pMotifRnd(idx,:);
for i = 1:3
    fprintf('%s %s %s\n',tbl.A{idx(i)},tbl.B{idx(i)},tbl.C{idx(i)});
end
fprintf('-----\n');
h = bar(1:15,[dev(:,1:15)]',1,'BaseValue',1,'EdgeColor','none'); 
h(1).FaceColor = 0.1.*ones(1,3);
h(2).FaceColor = 0.3.*ones(1,3);
h(3).FaceColor = 0.7.*ones(1,3);
set(gca,'TickDir','out','Box','off','XLim',[0.5 15.5],'XTick', ...
        [1:16],'YScale','log','Ylim',[min(dev(:)) max(dev(:))]); 
title('Layer 2/3');

% Columns
% C1-C2-C3 (19)
% D2-C2-E2 (20)
% C2-L4EXC C2-L4EXC C2-VPM (1794)
subplot(2,2,3);
idx = [19 20 1794];
dev = tbl.pMotif(idx,:)./tbl.pMotifRnd(idx,:);
for i = 1:3
    fprintf('%s %s %s\n',tbl.A{idx(i)},tbl.B{idx(i)},tbl.C{idx(i)});
end
fprintf('-----\n');
h = bar(1:15,[dev(:,1:15)]',1,'BaseValue',1,'EdgeColor','none'); 
h(1).FaceColor = 0.1.*ones(1,3);
h(2).FaceColor = 0.3.*ones(1,3);
h(3).FaceColor = 0.7.*ones(1,3);
set(gca,'TickDir','out','Box','off','XLim',[0.5 15.5],'XTick', ...
        [1:16],'YScale','log','Ylim',[min(dev(:)) max(dev(:))]); 
title('Columns');

set(f1,'PaperPositionMode', 'auto','Position',[0 0 600 300]);
print(f1,'-dsvg','-r600',[figurePath 'TripletMotifs_Examples_Supplements.svg']);