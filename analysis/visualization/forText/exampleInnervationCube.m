%% Numbers for Figure 2B,C,F 
% Output:
% Innervation: 0.0090
% 0: 99.11% (9.91e+01%)
% 1: 0.89% (8.87e-01%)
% 2: 0.00% (3.97e-03%)
% 3: 0.00% (1.19e-05%)
% Number of overlapping volumes = 13
% Max p across cubes = 6.21
% boutons in cube: 10.1814
% spines in cube: 113.6829
% all spines in cube: 129263.4346
% all exc spines in cube: 125370.8383
% all inh spines in cube: 3892.5963
% total spines L5PT: 23994.9457
% exc boutons(vS1) overlapping with L5PT: 33566426
% inh boutons(vS1) overlapping with L5PT: 5343557
% total boutons(vS1) overlapping with L5PT: 38909983
%
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
addpath([matlabPath 'functions\']);
dataPath1 = [matlabPath 'data\example_selected_cells\Innervation\'];
dataPath2 = [matlabPath 'data\amiraFiles\'];

%%
% example cube
% [-100 300 350] [-50 350 400]
% -> [-75 325 375]
targetCube = [-75 325 375];

[I, ~, ~, BB] = readAmira([dataPath1 'L2PY_L5PT_InnervationDensity.am']);
x = BB(1):50:BB(2);
y = BB(3):50:BB(4);
z = BB(5):50:BB(6);
idx = [find(x==targetCube(1)) find(y==targetCube(2)) ...
            find(z==targetCube(3))]; 
I_in_cube = I(idx(1),idx(2),idx(3));
numContacts = [0:1:3];
freqContacts = computeNumSynapses(I_in_cube,numContacts);

fprintf('Innervation: %.4f\n',I_in_cube); 
for i = 1:length(numContacts)
   fprintf('%d: %.2f%% (%.2e%%)\n',numContacts(i), ...
       freqContacts(i)*100,freqContacts(i)*100); 
end

fprintf('Number of overlapping volumes = %d\n',sum(I(:)>0));
fprintf('Max p across cubes = %.2f\n',max(1-exp(-I(:))).*100);

%%
% L5PT Spines in cube
[L5PT_spines, ~, ~, BB_spines] = readAmira([dataPath1 'L5PT_301854_spineDensity.am']);
x = BB_spines(1):50:BB_spines(2);
y = BB_spines(3):50:BB_spines(4);
z = BB_spines(5):50:BB_spines(6);
idx = [find(x==targetCube(1)) find(y==targetCube(2)) ...
            find(z==targetCube(3))]; 
L5PT_spine_in_cube = L5PT_spines(idx(1),idx(2),idx(3));

% L2PY Boutons in cube
[L2PY_boutons, ~, ~, BB] = readAmira([dataPath1 'L2PY_1012932_boutonDensity.am']);
x = BB(1):50:BB(2);
y = BB(3):50:BB(4);
z = BB(5):50:BB(6);
idx = [find(x==targetCube(1)) find(y==targetCube(2)) ...
            find(z==targetCube(3))]; 
L2PY_bouton_in_cube = L2PY_boutons(idx(1),idx(2),idx(3));

% All POSTs in cube
[all_spines, ~, ~, BB] = readAmira([dataPath2 'agg_pst_exc.am']);
x = BB(1):50:BB(2);
y = BB(3):50:BB(4);
z = BB(5):50:BB(6);
idx = [find(x==targetCube(1)) find(y==targetCube(2)) ...
            find(z==targetCube(3))]; 
all_spines_in_cube = all_spines(idx(1),idx(2),idx(3));

% All inh POSTs in cube
[tmp_spines, ~, ~, BB] = readAmira([dataPath2 'agg_pst_50-50-50_exc-to-inh.am']);
x = BB(1):50:BB(2);
y = BB(3):50:BB(4);
z = BB(5):50:BB(6);
idx = [find(x==targetCube(1)) find(y==targetCube(2)) ...
            find(z==targetCube(3))]; 
all_spines_in_cube_inh = tmp_spines(idx(1),idx(2),idx(3));

% All exc POSTs in cube
[tmp_spines, ~, ~, BB] = readAmira([dataPath2 'agg_pst_50-50-50_exc-to-exc.am']);
x = BB(1):50:BB(2);
y = BB(3):50:BB(4);
z = BB(5):50:BB(6);
idx = [find(x==targetCube(1)) find(y==targetCube(2)) ...
            find(z==targetCube(3))]; 
all_spines_in_cube_exc = tmp_spines(idx(1),idx(2),idx(3));

fprintf('boutons in cube: %.4f\n',L2PY_bouton_in_cube); 
fprintf('spines in cube: %.4f\n',L5PT_spine_in_cube); 
fprintf('all spines in cube: %.4f\n',all_spines_in_cube);
fprintf('all exc spines in cube: %.4f\n',all_spines_in_cube_exc); 
fprintf('all inh spines in cube: %.4f\n',all_spines_in_cube_inh); 
fprintf('total spines L5PT: %.4f\n',sum(L5PT_spines(:)));

%% All boutons that overlap with example cells
xSpines = BB_spines(1):50:BB_spines(2);
ySpines = BB_spines(3):50:BB_spines(4);
zSpines = BB_spines(5):50:BB_spines(6);

% EXC Boutons inside vS1 + VPM
[exc_boutons, ~, ~, BB] = readAmira([dataPath2 ...
                            'boutons_exc-inside_50-50-50_model-volume.am']);
                    
x = BB(1):50:BB(2);
y = BB(3):50:BB(4);
z = BB(5):50:BB(6);
numExcBoutonsOverlap = 0; 

for xx = 1:length(xSpines)
    for yy = 1:length(ySpines)
        for zz = 1:length(zSpines)
            if L5PT_spines(xx,yy,zz)>0
                xtmp = find(x==xSpines(xx));
                ytmp = find(y==ySpines(xx));
                ztmp = find(z==zSpines(xx));
                
                if isempty(xtmp) || isempty(ytmp) || isempty(ztmp)
                   warning('no bouton density found!'); 
                end
                
                numExcBoutonsOverlap = numExcBoutonsOverlap + ...
                                        exc_boutons(xtmp,ytmp,ztmp);
            end
        end
    end
end

% INH Boutons inside vS1
[inh_boutons, ~, ~, BB] = readAmira([dataPath2 ...
                        'boutons_inh-inside_50-50-50_model-volume.am']);
                    
x = BB(1):50:BB(2);
y = BB(3):50:BB(4);
z = BB(5):50:BB(6);
numInhBoutonsOverlap = 0; 

for xx = 1:length(xSpines)
    for yy = 1:length(ySpines)
        for zz = 1:length(zSpines)
            if L5PT_spines(xx,yy,zz)>0
                xtmp = find(x==xSpines(xx));
                ytmp = find(y==ySpines(xx));
                ztmp = find(z==zSpines(xx));
                
                if isempty(xtmp) || isempty(ytmp) || isempty(ztmp)
                   warning('no bouton density found!'); 
                end
                
                numInhBoutonsOverlap = numInhBoutonsOverlap + ...
                                        inh_boutons(xtmp,ytmp,ztmp);
            end
        end
    end
end


fprintf('exc boutons(vS1) overlapping with L5PT: %.0f\n',numExcBoutonsOverlap);
fprintf('inh boutons(vS1) overlapping with L5PT: %.0f\n',numInhBoutonsOverlap);
fprintf('total boutons(vS1) overlapping with L5PT: %.0f\n', ...
        numInhBoutonsOverlap+numExcBoutonsOverlap);