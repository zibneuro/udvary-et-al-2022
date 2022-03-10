%% Extract Synapse Density for comparison against measurements along D2 column
% Outputs: SynapseDensitiesAlongDepth.mat
% 
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';
inputPath = [matlabPath 'data\amiraFiles\'];
outputPath = [matlabPath 'data\'];
                       
cond = {'all','exc','inh'};
BarrelID = 'D2';
[BarrelCenter,BarrelzAxis,BarrelRadius] = getBarrelDefinition(BarrelID);
BarrelRadius = BarrelRadius + 50; 
barrelSeptumFactor = 1; 
D2offset = 707;
vol = 50*50*50; 
synDensity = struct; 
LayerBordersD2 = [0 157 296 575 900 1411 1973];

%%
for c = 1:3
    
    depth = []; 
    density = []; 
    
    % Load in Bouton Density
    filename = [inputPath 'boutons_' cond{c} ...
                    '-inside_50-50-50_model-volume.am'];
    [BoutonDensity, ~, ~, BB] = readAmira(filename);

    BBx = BB(1):50:BB(2);
    BBy = BB(3):50:BB(4);
    BBz = BB(5):50:BB(6);
    
    fid = fopen([outputPath cond{c} 'BoutonDensities.csv'],'w');
    fprintf(fid,'x,y,z,AxisDistance,Depth,#boutons,\n');
    
    for z = BBz
        
        depthtmp = D2offset - z; 
        
        % Only voxels that are in their entirty below pia
        % and in their entirty above the WM
        if ((depthtmp<LayerBordersD2(1)+50) || (depthtmp>LayerBordersD2(end)-50))
            continue;
        end
                               
        for x = BBx
            for  y = BBy
                
                dist = distPointToLine([x y z], ...
                            BarrelCenter,BarrelCenter+BarrelzAxis);

                if ((dist<=barrelSeptumFactor*BarrelRadius) && ...
                        BoutonDensity(x==BBx,y==BBy,z==BBz)>0)
                    
                    density = [density ...
                        BoutonDensity(x==BBx,y==BBy,z==BBz)/vol];
                    depth = [depth depthtmp]; 
                    
                    fprintf(fid,'%d,%d,%d,%d,%.2f,%.2f,\n',x,y,z, ...
                        dist,depthtmp,BoutonDensity(x==BBx,y==BBy,z==BBz));    
                end
            end
        end
    end
    
    eval(['synDensity.' cond{c} '_density = density;']);
    eval(['synDensity.' cond{c} '_depth = depth;']);
    fclose(fid);
end

save([outputPath 'SynapseDensitiesAlongDepth.mat'],'synDensity');