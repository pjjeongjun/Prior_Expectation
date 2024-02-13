%% Direction Tuning
% align preferred direction to 0 degree to combine data
% include only when there are tuning data of both two prior block
% Fig 4.
clear; clc; close all;
load('Tuning.mat');

%% fit tuning curve to each monkey's data

% choose a monkey
monkey = 'B'; % A for Monkey A, B for Monkey B
if strcmp(monkey,'A')
    errorBar = data.A.errorBar;
    mResponse = data.A.mResponse;
    stdResponse = data.A.stdResponse;
elseif strcmp(monkey,'B')
    errorBar = data.B.errorBar;
    mResponse = data.B.mResponse;
    stdResponse = data.B.stdResponse;
end

% convert tuning data from cell to array
mtuning1 = []; stuning1 = []; mtuning2 = []; stuning2 = []; mtuning3 = []; stuning3 = []; mtuning4 = []; stuning4 = [];
for i = 1:size(mResponse,2)

%     % For specific neurons based on the difference between preferred and prior dirs
%     load(['idx_' monkey '.mat']);
%     lidx for larger 30%, sidx for smaller 30%
%     for n = 1:length(lidx)
%         i = lidx(n); 

    mtmp1 = []; stmp1 = []; mtmp2 = []; stmp2 = []; mtmp3 = []; stmp3 = []; mtmp4 = []; stmp4 = [];
    for j = 1:12
        if ~isempty(mResponse{1,i})
            mtmp1 = [mtmp1; mResponse{1,i}{j}];
            stmp1 = [stmp1; stdResponse{1,i}{j}];
        end
        if ~isempty(mResponse{2,i})
            mtmp2 = [mtmp2; mResponse{2,i}{j}];
            stmp2 = [stmp2; stdResponse{2,i}{j}];
        end
        if ~isempty(mResponse{3,i})
            mtmp3 = [mtmp3; mResponse{3,i}{j}];
            stmp3 = [stmp3; stdResponse{3,i}{j}];
        end
        if ~isempty(mResponse{4,i})
            mtmp4 = [mtmp4; mResponse{4,i}{j}];
            stmp4 = [stmp4; stdResponse{4,i}{j}];
        end
    end
    mtuning1 = [mtuning1 mtmp1];
    mtuning2 = [mtuning2 mtmp2];
    mtuning3 = [mtuning3 mtmp3];
    mtuning4 = [mtuning4 mtmp4];
    
    stuning1 = [stuning1 stmp1];
    stuning2 = [stuning2 stmp2];
    stuning3 = [stuning3 stmp3];
    stuning4 = [stuning4 stmp4];
end
save(['' monkey '_data.mat'],'mtuning1','mtuning2','mtuning3','mtuning4','stuning1','stuning2','stuning3','stuning4');

mmtuning1 = mean(mtuning1,2);
mmtuning2 = mean(mtuning2,2);
mmtuning3 = mean(mtuning3,2);
mmtuning4 = mean(mtuning4,2);

mstuning1 = mean(stuning1,2);
mstuning2 = mean(stuning2,2);
mstuning3 = mean(stuning3,2);
mstuning4 = mean(stuning4,2);

%fitting
directions = [0 30 60 90 120 150 180 210 240 270 300 330]';
dirTuningInfo{1} = dirTuningEst(directions, mmtuning1, mstuning1); %b1c100
dirTuningInfo{2} = dirTuningEst(directions, mmtuning2, mstuning2); %b2c100
dirTuningInfo{3} = dirTuningEst(directions, mmtuning3, mstuning3); %b1c008
dirTuningInfo{4} = dirTuningEst(directions, mmtuning4, mstuning4); %b2c008

%draw figure

%gaussian
plot(dirTuningInfo{1}.estimates.gaussian.directions, dirTuningInfo{1}.estimates.gaussian.tuningCurve, 'k-.','color',[0 0 0]); %b1c100
hold on;
plot(dirTuningInfo{2}.estimates.gaussian.directions, dirTuningInfo{2}.estimates.gaussian.tuningCurve, 'k-.','color',[1 0 0]); %b2c100
hold on;
plot(dirTuningInfo{3}.estimates.gaussian.directions, dirTuningInfo{3}.estimates.gaussian.tuningCurve, 'k-.','color',[0.5 0.5 0.5]); %b1c008
hold on;
plot(dirTuningInfo{4}.estimates.gaussian.directions, dirTuningInfo{4}.estimates.gaussian.tuningCurve, 'k-.','color',[0.8 0.5 0.5]); %b2c008
hold on;

% %circular-gaussian
% plot(dirTuningInfo{1}.estimates.directions, dirTuningInfo{1}.estimates.tuningCurve, 'k-.','color',[0 0 0]); %b1c100
% hold on;
% plot(dirTuningInfo{2}.estimates.directions, dirTuningInfo{2}.estimates.tuningCurve, 'k-.','color',[1 0 0]); %b2c100
% hold on;
% plot(dirTuningInfo{3}.estimates.directions, dirTuningInfo{3}.estimates.tuningCurve, 'k-.','color',[0.5 0.5 0.5]); %b1c008
% hold on;
% plot(dirTuningInfo{4}.estimates.directions, dirTuningInfo{4}.estimates.tuningCurve, 'k-.','color',[0.8 0.5 0.5]); %b2c008
% hold on;

errorbar(directions, mmtuning1, mstuning1/sqrt(size(stuning1,2)), 'ro', 'markerfacecolor', [0 0 0], 'linewidth', 0.4, 'color', [0 0 0], 'markersize', 4);
hold on;
errorbar(directions, mmtuning2, mstuning2/sqrt(size(stuning2,2)), 'ro', 'markerfacecolor', [1 0 0], 'linewidth', 0.4, 'color', [1 0 0], 'markersize', 4);
hold on;
errorbar(directions, mmtuning3, mstuning3/sqrt(size(stuning3,2)), 'ro', 'markerfacecolor', [0.5 0.5 0.5], 'linewidth', 0.4, 'color', [0.5 0.5 0.5], 'markersize', 4);
hold on;
errorbar(directions, mmtuning4, mstuning4/sqrt(size(stuning4,2)), 'ro', 'markerfacecolor', [0.8 0.5 0.5], 'linewidth', 0.4, 'color', [0.8 0.5 0.5], 'markersize', 4);

savefig(['Tuning_' monkey '.fig']);

% t-test
ttestResult100 = []; ttestResult008 = [];
for k = 1:12
[h,p,ci,stats] = ttest(mtuning1(k,:),mtuning2(k,:));
ttestResult100(k,1) = h;
ttestResult100(k,2) = p;
ttestResult100(k,3) = stats.tstat;
ttestResult100(k,4) = nanmean(mtuning1(k,:));
ttestResult100(k,5) = nanmean(mtuning2(k,:));

[h,p,ci,stats] = ttest(mtuning3(k,:),mtuning4(k,:));
ttestResult008(k,1) = h;
ttestResult008(k,2) = p;
ttestResult008(k,3) = stats.tstat;
ttestResult008(k,4) = nanmean(mtuning3(k,:));
ttestResult008(k,5) = nanmean(mtuning4(k,:));
end

save(['stats_' monkey '.mat'],'ttestResult100','ttestResult008');

%% fit tuning curve to combined data
clear; clc;

% convert tuning data from cell to array
load('A_data.mat');
tmp1 = mtuning1; tmp2 = mtuning2; tmp3 = mtuning3; tmp4 = mtuning4; 
tmp5 = stuning1; tmp6 = stuning2; tmp7 = stuning3; tmp8 = stuning4;
load('B_data.mat');
mtuning1 = [tmp1 mtuning1]; mtuning2 = [tmp2 mtuning2]; mtuning3 = [tmp3 mtuning3]; mtuning4 = [tmp4 mtuning4];
stuning1 = [tmp5 stuning1]; stuning2 = [tmp6 stuning2]; stuning3 = [tmp7 stuning3]; stuning4 = [tmp8 stuning4];

mmtuning1 = mean(mtuning1,2);
mmtuning2 = mean(mtuning2,2);
mmtuning3 = mean(mtuning3,2);
mmtuning4 = mean(mtuning4,2);

mstuning1 = mean(stuning1,2);
mstuning2 = mean(stuning2,2);
mstuning3 = mean(stuning3,2);
mstuning4 = mean(stuning4,2);

%fitting
directions = [0 30 60 90 120 150 180 210 240 270 300 330]';
dirTuningInfo{1} = dirTuningEst(directions, mmtuning1, mstuning1); %b1c100
dirTuningInfo{2} = dirTuningEst(directions, mmtuning2, mstuning2); %b2c100
dirTuningInfo{3} = dirTuningEst(directions, mmtuning3, mstuning3); %b1c008
dirTuningInfo{4} = dirTuningEst(directions, mmtuning4, mstuning4); %b2c008

%draw figure

%gaussian
plot(dirTuningInfo{1}.estimates.gaussian.directions, dirTuningInfo{1}.estimates.gaussian.tuningCurve, 'k-.','color',[0 0 0]); %b1c100
hold on;
plot(dirTuningInfo{2}.estimates.gaussian.directions, dirTuningInfo{2}.estimates.gaussian.tuningCurve, 'k-.','color',[1 0 0]); %b2c100
hold on;
plot(dirTuningInfo{3}.estimates.gaussian.directions, dirTuningInfo{3}.estimates.gaussian.tuningCurve, 'k-.','color',[0.5 0.5 0.5]); %b1c008
hold on;
plot(dirTuningInfo{4}.estimates.gaussian.directions, dirTuningInfo{4}.estimates.gaussian.tuningCurve, 'k-.','color',[0.8 0.5 0.5]); %b2c008
hold on;

% %circular-gaussian
% plot(dirTuningInfo{1}.estimates.directions, dirTuningInfo{1}.estimates.tuningCurve, 'k-.','color',[0 0 0]); %b1c100
% hold on;
% plot(dirTuningInfo{2}.estimates.directions, dirTuningInfo{2}.estimates.tuningCurve, 'k-.','color',[1 0 0]); %b2c100
% hold on;
% plot(dirTuningInfo{3}.estimates.directions, dirTuningInfo{3}.estimates.tuningCurve, 'k-.','color',[0.5 0.5 0.5]); %b1c008
% hold on;
% plot(dirTuningInfo{4}.estimates.directions, dirTuningInfo{4}.estimates.tuningCurve, 'k-.','color',[0.8 0.5 0.5]); %b2c008
% hold on;

errorbar(directions, mmtuning1, mstuning1/sqrt(size(stuning1,2)), 'ro', 'markerfacecolor', [0 0 0], 'linewidth', 0.4, 'color', [0 0 0], 'markersize', 4);
hold on;
errorbar(directions, mmtuning2, mstuning2/sqrt(size(stuning2,2)), 'ro', 'markerfacecolor', [1 0 0], 'linewidth', 0.4, 'color', [1 0 0], 'markersize', 4);
hold on;
errorbar(directions, mmtuning3, mstuning3/sqrt(size(stuning3,2)), 'ro', 'markerfacecolor', [0.5 0.5 0.5], 'linewidth', 0.4, 'color', [0.5 0.5 0.5], 'markersize', 4);
hold on;
errorbar(directions, mmtuning4, mstuning4/sqrt(size(stuning4,2)), 'ro', 'markerfacecolor', [0.8 0.5 0.5], 'linewidth', 0.4, 'color', [0.8 0.5 0.5], 'markersize', 4);

savefig('Tuning_Monkey.fig');

% t-test
ttestResult100 = []; ttestResult008 = [];
for k = 1:12
    [h,p,ci,stats] = ttest(mtuning1(k,:),mtuning2(k,:));
    ttestResult100(k,1) = h;
    ttestResult100(k,2) = p;
    ttestResult100(k,3) = stats.tstat;
    ttestResult100(k,4) = nanmean(mtuning1(k,:));
    ttestResult100(k,5) = nanmean(mtuning2(k,:));

    [h,p,ci,stats] = ttest(mtuning3(k,:),mtuning4(k,:));
    ttestResult008(k,1) = h;
    ttestResult008(k,2) = p;
    ttestResult008(k,3) = stats.tstat;
    ttestResult008(k,4) = nanmean(mtuning3(k,:));
    ttestResult008(k,5) = nanmean(mtuning4(k,:));
end

save('stats_Monkey.mat','ttestResult100','ttestResult008');
