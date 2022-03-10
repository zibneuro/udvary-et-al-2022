%% Figure S3B
% Author: Daniel Udvary (Max Planck Institute for Neurobiology of Behavior – caesar)
clear all
close all
clc

matlabPath = 'D:\udvary-et-al-2022\analysis\';

%%
L5PT_ID = 301854;
L2PY_ID = 748854;
L6ACC_ID = 199678;

%% L2PY->L5PT
filename = [matlabPath 'data\CellMatrix_C2\L2PY_L5PT.mat'];
load(filename); 
idxPre = (I.PreCellID==L2PY_ID);
idxPost = (I.PostCellID==L5PT_ID);
Ival = I.I(idxPre,idxPost);
fprintf('L2PY->L5PT: %.4f (p = %.4f%%)\n',Ival,100*(1-exp(-Ival)));

%% L5PT->L2PY
filename = [matlabPath 'data\CellMatrix_C2\L5PT_L2PY.mat'];
load(filename); 
idxPre = (I.PreCellID==L5PT_ID);
idxPost = (I.PostCellID==L2PY_ID);
Ival = I.I(idxPre,idxPost);
fprintf('L5PT->L2PY: %.4f (p = %.4f%%)\n',Ival,100*(1-exp(-Ival)));

%% L2PY->L6ACC
filename = [matlabPath 'data\CellMatrix_C2\L2PY_L6ACC.mat'];
load(filename); 
idxPre = (I.PreCellID==L2PY_ID);
idxPost = (I.PostCellID==L6ACC_ID);
Ival = I.I(idxPre,idxPost);
fprintf('L2PY->L6ACC: %.4f (p = %.4f%%)\n',Ival,100*(1-exp(-Ival)));

%% L6ACC->L2PY
filename = [matlabPath 'data\CellMatrix_C2\L6ACC_L2PY.mat'];
load(filename); 
idxPre = (I.PreCellID==L6ACC_ID);
idxPost = (I.PostCellID==L2PY_ID);
Ival = I.I(idxPre,idxPost);
fprintf('L6ACC->L2PY: %.4f (p = %.4f%%)\n',Ival,100*(1-exp(-Ival)));

%% L5PT->L6ACC
filename = [matlabPath 'data\CellMatrix_C2\L5PT_L6ACC.mat'];
load(filename); 
idxPre = (I.PreCellID==L5PT_ID);
idxPost = (I.PostCellID==L6ACC_ID);
Ival = I.I(idxPre,idxPost);
fprintf('L5PT->L6ACC: %.4f (p = %.4f%%)\n',Ival,100*(1-exp(-Ival)));

%% L6ACC->L5PT
filename = [matlabPath 'data\CellMatrix_C2\L6ACC_L5PT.mat'];
load(filename); 
idxPre = (I.PreCellID==L6ACC_ID);
idxPost = (I.PostCellID==L5PT_ID);
Ival = I.I(idxPre,idxPost);
fprintf('L6ACC->L5PT: %.4f (p = %.4f%%)\n',Ival,100*(1-exp(-Ival)));