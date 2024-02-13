% This code is a modified version of the original code from Mante et al 2013.
% V. Mante, D. Sussillo, K. V. Shenoy, W. T. Newsome, 
% Context-dependent computation by recurrent dynamics in prefrontal cortex. Nature. 503, 78–84 (2013).

%% Add paths

clear; clc;

% TDR path
tdrDir = ''; % put your code directory
addpath(fullfile(tdrDir,'TDR_from_Mante_et_al_2013'));
addpath(fullfile(tdrDir,'TDR_from_Mante_et_al_2013','nansuite'));
addpath(fullfile(tdrDir,'TDR_from_Mante_et_al_2013','tools'));

%% Load data
% firing rate
% time = 1:20:981 (400 = stimulus onset)
% time window = 20ms(0~19ms) from time point

animal = 'A'; % A or B
load(['data_' animal '.mat']);

%% Parameters

plotflag = 1;

%% Condition averaged responses

% The conditions to use
task_index = [];
task_index.stim_dir = [1 1 2 2 3 3 3 3 4 4 5 5]'; %-120,-15,0,15,120:1,2,3,4,5
task_index.prior =    [2 2 1 1 2 2 1 1 1 1 2 2]'; %narrow:1,wide:2
task_index.contrast = [1 2 1 2 1 2 1 2 1 2 1 2]'; %low:1,high:2

% Condition average
dataC = tdrAverageCondition(dataT,task_index);

%% Process condition averaged responses

% Averaging parameters
avgpars = [];
avgpars.trial = [];
avgpars.time = [];

% Mean and STD across time and conditions
[meanC,stdC] = tdrMeanAndStd(dataC,avgpars);

% Normalization parameters
nrmlpars = [];
nrmlpars.ravg = meanC;
nrmlpars.rstd = stdC;
% nrmlpars.cnst = median(stdC)/2; % arbitrary

% Normalize
dataC_nrml = tdrNormalize(dataC,nrmlpars);

%% Linear regression

% Averaging parameters
avgpars = [];
avgpars.trial = [];
avgpars.time = [];

% Mean and STD across time and trials
[meanT,stdT] = tdrMeanAndStd(dataT,avgpars);

% Normalization parameters
nrmlpars = [];
nrmlpars.ravg = meanT;
nrmlpars.rstd = stdT;
% nrmlpars.cnst = median(stdT)/2; % arbitrary

% Normalize
dataT_nrml = tdrNormalize(dataT,nrmlpars);

% Regression parameters
regpars = [];
regpars.regressor = {...
    'b0';...
    'stim_dir';'prior';'contrast'};
regpars.regressor_normalization = 'max_abs';

% Linear regression
coef_fulUN = tdrRegression(dataT_nrml,regpars,plotflag);


%% Principal component analysis

% PCA parameters
pcapars = [];
% pcapars.trial_pca = dataC_nrml.task_index.diff==1;
pcapars.trial_prj = [];
pcapars.time_pca = [];
pcapars.plot_dimensions = 1:20;

% Compute PCA
[data_fulPC,fulUN_fulPC,varPC] = tdrPca(dataC_nrml,pcapars,plotflag);


%% Define mid-dimensional subspace

% Principle components to keep
pc2keep = 1:12;

% Projection matrix full space (unit basis) into PC subspace (PC basis)
midPC_fulUN = fulUN_fulPC(:,pc2keep)';


%% Variance explained by PCs

% Variance parameters
varpars = [];
varpars.time_var = [];
varpars.dim_sub = pc2keep;

% Compute variance
var_mid = tdrSubspaceVariance(data_fulPC,dataC_nrml,varpars,plotflag);


%% Denoise and smooth regression coefficients

% Subspace parameters
subpars = [];
subpars.subSUB_fulORG = midPC_fulUN;
subpars.dimension = data_fulPC.dimension(pc2keep);

% Project coefficients into subspace
[coef_midPC,coef_midUN] = tdrSubspaceProjection(coef_fulUN,subpars);

% Smoothing parameters
smthpars = [];
smthpars.filter = 'gauss'; % {'gauss';'box'}
smthpars.width = 0.04;

% Smooth coefficients
coef_midPC = tdrTemporalSmoothing(coef_midPC,smthpars);
coef_midUN = tdrTemporalSmoothing(coef_midUN,smthpars);


%% Temporal dynamics of regression coefficients

% Correlogram parameters
plotpars = [];
plotpars.name = {'stim_dir';'prior';'contrast'};
plotpars.plotpairs = 1;

% Coefficient correlogram
[~,~,h] = tdrVectorDynamics(coef_midUN,plotpars,plotflag);


%% Define regression vectors

% Regression vector parameters
vecpars = [];
vecpars.stim_dir.time_win = [];
vecpars.prior.time_win    = [];
vecpars.contrast.time_win = [];

% Compute regression vectors
vBeta = tdrVectorTimeAverage(coef_midUN,vecpars,plotflag);


%% Define task-related axes (orthogonalize regression vectors)

% Regression axes parameters
ortpars = [];
ortpars.name = {'stim_dir';'prior';'contrast'};

% Compute regression axes
[vAxes,lowUN_lowTA] = tdrVectorOrthogonalize(vBeta,ortpars);

% Projection matrix full space (unit basis) into task subspace (task basis)
lowTA_lowUN = lowUN_lowTA';


%% Responses in task-related subspace

% Subspace parameters
subpars = [];
subpars.subSUB_fulORG = lowTA_lowUN;
subpars.dimension = vAxes.name;

% Project responses into subspace
[dataC_lowTA,dataC_lowUN,varTA] = tdrSubspaceProjection(dataC_nrml,subpars,0);

% Variance parameters
varpars = [];
varpars.time_var = pcapars.time_pca;
var_low = tdrSubspaceVariance(dataC_lowTA,dataC_nrml,varpars,plotflag);


%% Plot 2D trajectories - condition averaged responses

% PICK AXES
plotflag = 1; % 1 = draw figure / 0 = no figure

%--- Data to plot ---
data = dataC_lowTA;

if plotflag
    
    % Common plot parameters
    % Markersize and linewidth. Use default if empty.
    plotpars = [];
    plotpars.markersize = [];
    plotpars.linewidth = [];
    plotpars.handle.figure = [];
    plotpars.handle.axes = [];
    plotpars.handle.plot = [];
    plotpars.average_trials = 0; % 1 or 0
    plotpars.type = '-o';
    plotpars.dataaspectratio = [1 1 1];
    plotpars.plotboxaspectratio = [1 1 1];
    plotpars.jtime_plot = 1:20; % time window (rf: data.time)
    plotpars.bootstrapN = [];

    %%%%%%%% dimension %%%%%%%%
    plotpars.dimension = {'prior','stim_dir'};

    %--- High contrast ---%
    plotpars.handle = [];
    plotpars.task_index.contrast = 2; %high
    plotpars.title = 'High contrast';

    % Choice 1
    plotpars.colormap = 'gray';
    plotpars.task_index.prior    = 2; %wide
    plotpars.task_index.stim_dir = [1 3 5];
    plotpars.markerface = 'full';
    [sdata_h_w, plotpars.handle] = tdrPlotResponse2D(data,plotpars);
    
    % Choice 2
    plotpars.colormap = 'blue';
    plotpars.task_index.prior = 1;    %narrow
    plotpars.task_index.stim_dir = [2 3 4];
    plotpars.legend = {'-120','0','120','-15','0','15'};
    plotpars.markerface = 'empty';
    [sdata_h_n, plotpars.handle] = tdrPlotResponse2D(data,plotpars);

    %--- Low contrast ---%
    plotpars.handle = [];
    plotpars.task_index.contrast = 1; %low
    plotpars.title = 'Low contrast';
    
    % Choice 1
    plotpars.colormap = 'gray';
    plotpars.task_index.prior = 2;    %wide
    plotpars.task_index.stim_dir = [1 3 5];
    plotpars.markerface = 'full';
    [sdata_l_w, plotpars.handle] = tdrPlotResponse2D(data,plotpars);
    
    % Choice 2
    plotpars.colormap = 'blue';
    plotpars.task_index.prior = 1;    %narrow
    plotpars.task_index.stim_dir = [2 3 4];
    plotpars.legend = {'-120','0','120','-15','0','15'};
    plotpars.markerface = 'empty';
    [sdata_l_n, plotpars.handle] = tdrPlotResponse2D(data,plotpars);
end
