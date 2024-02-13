%% Add paths

clear; clc;

% TDR path
tdrDir = '/Users/jeongjunpark/Desktop/Joonyeol/Mante';
addpath(fullfile(tdrDir,'TDR'));
addpath(fullfile(tdrDir,'TDR','nansuite'));
addpath(fullfile(tdrDir,'TDR','tools'));

%% Load data
% firing rate
% time = 101:20:621 (300 = stimulus onset)
% time window = 20ms(0~19ms) from time point

animal = 'A'; % A or B
load(['data_' animal '.mat']);

%% Parameters

plotflag = 1;

%% Condition averaged responses

% The conditions to use
task_index = [];
task_index.stim_dir = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 7 8 8 8 8 9 9 9 9 10 10 10 10 11 11 11 11 12 12 12 12]'; %0,30,...,300,330:1~12
task_index.prior =    [1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2 1 1 2 2]'; %narrow:1,wide:2
task_index.contrast = [1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2]'; %low:1,high:2

% Condition average
dataC = tdrAverageCondition(dataT,task_index);

% Remove cells including nan values
idx = [];
for n = 1:size(dataC.response,1)
    if ~isempty(find(isnan(dataC.response(n,:,:))))
        idx = [idx; n];
    end
end
dataC.response(idx,:,:)=[];
dataC.n_trial(idx,:)=[];
dataC.dimension(idx)=[];
dataT.unit(idx)=[];

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
    'prior';'contrast'};
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
plotpars.name = {'prior';'contrast'};
plotpars.plotpairs = 1;

% Coefficient correlogram
[~,~,h] = tdrVectorDynamics(coef_midUN,plotpars,plotflag);


%% Define regression vectors

% Regression vector parameters
vecpars = [];
vecpars.prior.time_win    = [];
vecpars.contrast.time_win = [];

% Compute regression vectors
vBeta = tdrVectorTimeAverage(coef_midUN,vecpars,plotflag);


%% Define task-related axes (orthogonalize regression vectors)

% Regression axes parameters
ortpars = [];
ortpars.name = {'prior';'contrast'};

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
%     plotpars.jtime_plot = 6:10; % time window (rf: data.time)
    plotpars.bootstrapN = [];

%     %%%%%%%% dimension %%%%%%%%
%     plotpars.dimension = {'prior','contrast'};
% 
%     %--- High contrast ---%
%     plotpars.handle = [];
%     plotpars.task_index.contrast = 2; %high
%     plotpars.title = 'High contrast';
% 
%     % Choice 1
%     plotpars.colormap = 'gray';
%     plotpars.task_index.prior    = 2; %wide
%     plotpars.task_index.stim_dir = 1:12;
%     plotpars.markerface = 'full';
%     [sdata_h_w, plotpars.handle] = tdrPlotResponse2D(data,plotpars);
%     
%     % Choice 2
%     plotpars.colormap = 'blue';
%     plotpars.task_index.prior = 1;    %narrow
%     plotpars.task_index.stim_dir = 1:12;
%     plotpars.legend = {'0','30','60','90','120','150','180','210','240','270','300','330','0','30','60','90','120','150','180','210','240','270','300','330'};
%     plotpars.markerface = 'empty';
%     [sdata_h_n, plotpars.handle] = tdrPlotResponse2D(data,plotpars);
% 
%     %--- Low contrast ---%
%     plotpars.handle = [];
%     plotpars.task_index.contrast = 1; %low
%     plotpars.title = 'Low contrast';
%     
%     % Choice 1
%     plotpars.colormap = 'gray';
%     plotpars.task_index.prior = 2;    %wide
%     plotpars.task_index.stim_dir = 1:12;
%     plotpars.markerface = 'full';
%     [sdata_l_w, plotpars.handle] = tdrPlotResponse2D(data,plotpars);
%     
%     % Choice 2
%     plotpars.colormap = 'blue';
%     plotpars.task_index.prior = 1;    %narrow
%     plotpars.task_index.stim_dir = 1:12;
%     plotpars.legend = {'0','30','60','90','120','150','180','210','240','270','300','330','0','30','60','90','120','150','180','210','240','270','300','330'};
%     plotpars.markerface = 'empty';
%     [sdata_l_n, plotpars.handle] = tdrPlotResponse2D(data,plotpars);

    %%%%%%%% dimension %%%%%%%%
    for n = 1:12
        if n == 1
            plotpars.dimension = {'prior','contrast'};
            plotpars.handle = [];
            plotpars.title = 'subspace';
        end

        plotpars.task_index.stim_dir = n;

        % Choice 1
        plotpars.task_index.contrast = 2; %high
        plotpars.colormap = 'blue';
        plotpars.task_index.prior    = [1 2];%narrow,wide
        plotpars.markerface = 'full';
        plotpars.dname = ['' animal '_high'];
        [sdata_h, plotpars.handle] = tdrPlotResponse2D(data,plotpars);

        % Choice 2
        plotpars.task_index.contrast = 1; %low
        plotpars.colormap = 'darkgray';
        plotpars.task_index.prior = [1 2];
        plotpars.dname = ['' animal '_low'];
        [sdata_l, plotpars.handle] = tdrPlotResponse2D(data,plotpars);
       
        sdata_h_all{n} = sdata_h;
        sdata_l_all{n} = sdata_l;
    end
end

% Concatenate all the dirs & Save x,y plot data
sdata = struct('high',[],'low',[]);
sdata.high = struct('wide',[],'narrow',[]);
sdata.high.wide = struct('x',[],'y',[]);
sdata.high.narrow = struct('x',[],'y',[]);
sdata.low = struct('wide',[],'narrow',[]);
sdata.low.wide = struct('x',[],'y',[]);
sdata.low.narrow = struct('x',[],'y',[]);

for n = 1:12
    sdata.high.wide.x = [sdata.high.wide.x; sdata_h_all{n}{2,1}];
    sdata.high.wide.y = [sdata.high.wide.y; sdata_h_all{n}{2,2}];

    sdata.high.narrow.x = [sdata.high.narrow.x; sdata_h_all{n}{1,1}];
    sdata.high.narrow.y = [sdata.high.narrow.y; sdata_h_all{n}{1,2}];

    sdata.low.wide.x = [sdata.low.wide.x; sdata_l_all{n}{2,1}];
    sdata.low.wide.y = [sdata.low.wide.y; sdata_l_all{n}{2,2}];

    sdata.low.narrow.x = [sdata.low.narrow.x; sdata_l_all{n}{1,1}];
    sdata.low.narrow.y = [sdata.low.narrow.y; sdata_l_all{n}{1,2}];
end
save(['subspace_position_' animal '.mat'],'sdata');