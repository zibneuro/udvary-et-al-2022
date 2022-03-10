%% Figure S2B
% Pathlength within vS1 model
%   by exc/inh/vpm and dend/axon
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
figuresPath = [matlabPath 'output\Figures\']; 
tablePath = [matlabPath 'output\Tables\'];

%%
str = fileread([matlabPath 'data\length_2021-05-04.json']);
data = jsondecode(str);

% Stack
% Length
lenExcAxon = str2double(data.length_axon_exc); % in um
lenInhAxon = str2double(data.length_axon_inh);
lenVPMAxon = str2double(data.length_VPM);
lenExcDend = str2double(data.length_dendrite_exc);
lenInhDend = str2double(data.length_dendrite_inh);

% Boutons
boutonsVPM = str2double(data.boutons_VPM);
boutonsInh = str2double(data.boutons_inh);
boutonsExc = str2double(data.boutons_exc);

% Post when presynapse is exc.
postExcToExc = str2double(data.pst_exc_to_exc);
postExcToInh = str2double(data.pst_exc_to_inh);

% Post when presynapse is inh.
postInhToExc = str2double(data.pst_inh_to_exc);
postInhToInh = str2double(data.pst_inh_to_inh);

% Length stats
len_label = {'ExcAxonLen','InhAxonLen','VPMAxonLen','ExcDendLen','InhDendLen'};
len_X = [lenExcAxon-lenVPMAxon lenInhAxon lenVPMAxon lenExcDend lenInhDend];
len_X = len_X./1e9; % in km

for i = 1:length(len_X)
   fprintf('%s: %.3fkm\n',len_label{i},len_X(i)); 
end
fprintf('TotalLen: %.3fkm (Axon: %.3fkm Dend: %.3fkm)\n',sum(len_X), ...
    (lenExcAxon+lenInhAxon)/1e9,(lenExcDend+lenInhDend)/1e9); 

% per mm³
vol_vS1 = 6.7; % 6.7 according to Amira surface
% cubes considered 6.65
fprintf('TotalLen: %.3fkm (Axon: %.3fkm Dend: %.3fkm) per mm³\n',sum(len_X)/vol_vS1, ...
    (lenExcAxon+lenInhAxon)/1e9/vol_vS1,(lenExcDend+lenInhDend)/1e9/vol_vS1); 
fprintf('-----------\n');
% Bouton stats
boutons_label = {'ExcBoutons','InhBoutons','VPMBoutons'};
boutons_X = [boutonsExc-boutonsVPM boutonsInh boutonsVPM];

for i = 1:length(boutons_X)
   fprintf('%s: %.0f\n',boutons_label{i},boutons_X(i)); 
end
fprintf('TotalBoutons: %.0f\n',sum(boutons_X)); 
fprintf('-----------\n');

tbl = table();
tbl.feature = {'path_length_excitatory_axons[km]', ...
                'path_length_inhibitory_axons[km]', ...
                'path_length_VPM_axons[km]', ...
                'path_length_excitatory_dendrites[km]',...
                'path_length_inhibitory_axons[km]', ...
                'excitatory_boutons','inhibitory_boutons','VPM_boutons'}';
tbl.value = [len_X boutons_X]';
writetable(tbl,[tablePath 'stats_of_barrel_cortex_structure.csv']);

% Postsynaptic sites stats when presynapse is excitatory
fprintf('PostExcToExc: %.0f\n',postExcToExc); 
fprintf('PostExcToInh: %.0f\n',postExcToInh); 
fprintf('TotalExcPostTargets: %.0f\n',postExcToExc+postExcToInh); 
fprintf('Discrepancy: ExcBoutons/SpineTargets: %.4f\n', ...
            boutonsExc/(postExcToExc+postExcToInh));
fprintf('-----------\n');
  
% Postsynaptic sites stats when presynapse is inhibitory
fprintf('PostInhToExc: %.0f\n',postInhToExc); 
fprintf('PostInhToInh: %.0f\n',postInhToInh); 
fprintf('TotalInhPostTargets: %.0f\n',postInhToExc+postInhToInh); 
fprintf('-----------\n');
fprintf('Discrepancy: ExcBoutons/SpineTargets: %.4f\n', ...
            boutonsInh/(postInhToExc+postInhToInh));      


%%
f1 = figure(1);
clf;
subplot(3,1,1);
b1 = barh(1,len_X,'stacked');
b1(1).FaceColor = [0.5 0.5 0.5]; % EXC Axon
b1(2).FaceColor = [1 1 1]; % INH Axon
b1(3).FaceColor = [0 0 0]; % VPM Axon
b1(4).FaceColor = [0.5 0.5 0.5]; % EXC Dend
b1(5).FaceColor = [1 1 1]; % INH Dend
set(gca,'Box','off','TickDir','out','YTick',[])
xlabel('path length (km)');

subplot(3,1,2);
b1 = barh(1,boutons_X./1e9,'stacked');
b1(1).FaceColor = [0.5 0.5 0.5]; % EXC boutons
b1(2).FaceColor = [1 1 1]; % INH boutons
b1(3).FaceColor = [0 0 0]; % VPM boutons
set(gca,'Box','off','TickDir','out','YTick',[])
xlabel('boutons (x10^9)');

subplot(3,1,3);
b1 = barh(1,[postExcToExc postExcToInh]./1e9,'stacked');
b1(1).FaceColor = [0.5 0.5 0.5]; % EXC Posts
b1(2).FaceColor = [1 1 1]; % INH Posts
set(gca,'Box','off','TickDir','out','YTick',[])
xlabel('posts (x10^9)');

set(f1,'PaperPositionMode','auto','Position',[0 0 300 300]); 
print(f1,'-painters','-dsvg','-r600',[figuresPath 'PathLength_vS1.svg']); 