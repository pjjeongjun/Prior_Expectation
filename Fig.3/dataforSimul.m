% Experimental data for simulation
% Fig. 3
%
% go to 'SVM.m' or 'PVD_and_MLE.m' to simulate and decode data

%% Data of one monkey using Circular Gaussian fitting
clear; clc; close all;

load('neural.mat');

%% Reshape data of each monkey

for m = 1:2
    Data = [];
    
    if m == 1
        monkey = 'A';
        mRespPfHigh = mRespPfHigh_A; mRespNullHigh = mRespNullHigh_A; mRespPfLow = mRespPfLow_A; mRespNullLow = mRespNullLow_A; sigmac100 = sigmac100_A; sigmac008 = sigmac008_A; cparams008 = cparams008_A; cparams100 = cparams100_A;
    else
        monkey = 'B';
        mRespPfHigh = mRespPfHigh_B; mRespNullHigh = mRespNullHigh_B; mRespPfLow = mRespPfLow_B; mRespNullLow = mRespNullLow_B; sigmac100 = sigmac100_B; sigmac008 = sigmac008_B; cparams008 = cparams008_B; cparams100 = cparams100_B;
    end
    
    % 1. Two contrasts (prior blocks are combined)
    Data(1).Contrast = 'High';
    Data(2).Contrast = 'Low';
    
    % 2. Peak and Null resp    
    mRespTunfittedH = []; tmp1 = []; tmp2 = []; % High
    for g = 1:size(mRespPfHigh,2)
        tmp1 = [tmp1 mRespPfHigh{g}];
        tmp2 = [tmp2 mRespNullHigh{g}];
    end
    mRespTunfittedH = [tmp1; tmp2];
    
    mRespTunfittedL = []; tmp1 = []; tmp2 = []; % Low
    for g = 1:size(mRespPfLow,2)
        tmp1 = [tmp1 mRespPfLow{g}];
        tmp2 = [tmp2 mRespNullLow{g}];
    end
    mRespTunfittedL = [tmp1; tmp2];
    
    Data(1).Tuning_pfnull = mRespTunfittedH;
    Data(2).Tuning_pfnull = mRespTunfittedL;
    
    % change num to nan (if prefer is smaller than null)
    idxH = []; % High
    for k = 1:size(Data(1).Tuning_pfnull,2)
        if ~isnan(Data(1).Tuning_pfnull(1,k))
            if Data(1).Tuning_pfnull(1,k)-Data(1).Tuning_pfnull(2,k) < 0
                Data(1).Tuning_pfnull(1,k) = nan;
                Data(1).Tuning_pfnull(2,k) = nan;
                idxH = [idxH k];
            end
        end
    end
    
    idxL = []; %Low
    for k = 1:size(Data(2).Tuning_pfnull,2)
        if ~isnan(Data(2).Tuning_pfnull(1,k))
            if Data(2).Tuning_pfnull(1,k)-Data(2).Tuning_pfnull(2,k) < 0
                Data(2).Tuning_pfnull(1,k) = nan;
                Data(2).Tuning_pfnull(2,k) = nan;
                idxL = [idxL k];
            end
        end
    end
    
    % 3. Sigma for tuning width
    sig100 = []; sig008 = [];
    for s = 1:size(sigmac100,2)
        if ~isempty(sigmac100{s})
            sig100 = [sig100; sigmac100{s}];
        else
            sig100 = [sig100; nan];
        end
        if ~isempty(sigmac008{s})
            sig008 = [sig008; sigmac008{s}];
        else
            sig008 = [sig008; nan];
        end
    end
    Data(1).Sigma = sig100;
    Data(2).Sigma = sig008;

    % 4. Fano factor     
    if strcmp(monkey,'A')
        Data(1).Fanofactor = targetHigh_A; % High
        Data(2).Fanofactor = targetLow_A; % Low
    elseif strcmp(monkey,'B')
        Data(1).Fanofactor = targetHigh_B; % High
        Data(2).Fanofactor = targetLow_B; % Low
    end
    
    % 5. Neuron-Neuron correlation
    if strcmp(monkey,'A')
        Data(1).NNcorr = nnCorr_high_a; %High
        Data(2).NNcorr = nnCorr_low_a; %Low
    elseif strcmp(monkey,'B')
        Data(1).NNcorr = nnCorr_high_b; %High
        Data(2).NNcorr = nnCorr_low_b; %Low
    end
    
    % 6. Circular Gaussian parameters
    tmp100 = [];
    for i = 1:length(cparams100)
        if isempty(cparams100{i})
            tmp1 = [NaN NaN NaN];
            tmp100 = [tmp100; tmp1];
        else
            tmp2 = [cparams100{i}(1) cparams100{i}(2) cparams100{i}(3)];
            tmp100 = [tmp100; tmp2];
        end
    end
    Data(1).cparams = tmp100;
    
    tmp008 = [];    
    for i = 1:length(cparams008)
        if isempty(cparams008{i})
            tmp1 = [NaN NaN NaN];
            tmp008 = [tmp008; tmp1];
        else
            tmp2 = [cparams008{i}(1) cparams008{i}(2) cparams008{i}(3)];
            tmp008 = [tmp008; tmp2];
        end
    end
    Data(2).cparams = tmp008;
    
    save(['' monkey '.mat'],'Data');
end

%% combine two monkeys data

load('A.mat');
Ace = Data;
load('B.mat');
Bingo = Data;

Data = [];

% contrast type name
Data(1).Contrast = 'High';
Data(2).Contrast = 'Low';

% peak and null
Data(1).Tuning_pfnull = [Ace(1).Tuning_pfnull Bingo(1).Tuning_pfnull];
Data(2).Tuning_pfnull = [Ace(2).Tuning_pfnull Bingo(2).Tuning_pfnull];

% sigma
Data(1).Sigma = [Ace(1).Sigma; Bingo(1).Sigma];
Data(2).Sigma = [Ace(2).Sigma; Bingo(2).Sigma];

% Fano factor
Data(1).Fanofactor = [Ace(1).Fanofactor; Bingo(1).Fanofactor];
Data(2).Fanofactor = [Ace(2).Fanofactor; Bingo(2).Fanofactor];

% Neuron-Neuron correlation
Data(1).NNcorr = [Ace(1).NNcorr; Bingo(1).NNcorr];
Data(2).NNcorr = [Ace(2).NNcorr; Bingo(2).NNcorr];

% regression
% wStep = [0:5:200]; 18 = 85ms.
Data(1).Regression = rvH;
Data(2).Regression = rvL;

% circular gaussian params
Data(1).cparams = [Ace(1).cparams; Bingo(1).cparams];
Data(2).cparams = [Ace(2).cparams; Bingo(2).cparams];

save('Monkeys_cir.mat','Data');
