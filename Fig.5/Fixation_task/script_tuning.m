% Spontaneous neural activity during fixation task
% Fig. 5d, fig. S7d, and fig. S9, b,d,f 

%% PSTH of spontaneous activity
clear; clc; close all;

evCrit = 0.5; % explained variance of fitted tuning (> 50%)
monkey = 'A'; % A for Monkey A, B for Monkey B

load(['Data_' monkey '_fixation_ev' num2str(evCrit) '.mat']);
wide = []; narrow = [];
for i = 1:12
    wide = cat(3,wide,b1c100{i});
    wide = cat(3,wide,b1c008{i});
    narrow = cat(3,narrow,b2c100{i});
    narrow = cat(3,narrow,b2c008{i});
end

wide = nanmean(wide,3);
narrow = nanmean(narrow,3);

% Gaussian smoothing
for i = 1:size(wide,1)
    wide(i,:) = gaussianFilter(wide(i,:),10);
end
for i = 1:size(narrow,1)
    narrow(i,:) = gaussianFilter(narrow(i,:),10);
end

% PSTH
figure;
x = -299:0;
errorbar(x,nanmean(wide),nanstd(wide)/sqrt(size(wide,1)));
hold on;
x = -299:0;
errorbar(x,nanmean(narrow), nanstd(narrow)/sqrt(size(narrow,1)));
legend('wide-prior','narrow-prior');
xlabel('Time (ms)'); ylabel('Firing rate'); title('Fixation task');
savefig(['psth_' monkey '.fig'])
save(['' monkey '.mat'],'wide','narrow');

% Mean spontaneous activity
wide = nanmean(wide(:,200:end),2);
narrow = nanmean(narrow(:,200:end),2);

% ttest
ttestResult = [];
[h,p,ci,stats] = ttest(wide,narrow);
ttestResult(1,1) = h;
ttestResult(1,2) = p;
ttestResult(1,3) = stats.tstat;
ttestResult(1,4) = nanmean(wide);
ttestResult(1,5) = nanmean(narrow);

fig = figure;
c = categorical({'wide-prior','narrow-prior'});
m = [nanmean(wide),nanmean(narrow)];
sd = [nanstd(wide)/sqrt(length(~isnan(wide))),nanstd(narrow)/sqrt(length(~isnan(narrow)))];
errorbar(c,m,sd);
hold on;
bar(c(1),m(1),'b')
hold on;
bar(c(2),m(2),'r')
title(['Mean spontaneous activity, p = ' num2str(ttestResult(1,2)) '']); ylabel('Firing rate');
savefig(['bar_' monkey '.fig'])
save(['' monkey '_mean.mat'],'wide','narrow');

%% Mean spontaneous activity - delta tuning
clear; clc;

load('delta_theta.mat');
load('A_mean.mat');
wide_a = wide; narrow_a = narrow;
load('B_mean.mat');
wide_b = wide; narrow_b = narrow;

%concatenate
wide = [wide_a; wide_b];
narrow = [narrow_a; narrow_b];
dirD = [dirD_a; dirD_b];

%regression - log ratio
[r,int] = regress(log(narrow./wide),[ones(size(dirD,1),1) dirD]);
y = r(2)*dirD+r(1);

%correlation
idx = find(~isnan(dirD)&~isnan(wide)&~isnan(narrow));
c1 = log(narrow(idx)./wide(idx));
c2 = dirD(idx);
[rho,p] = corr(c1,c2,'type','Spearman');

figure;
scatter(dirD_a, log(narrow_a./wide_a),'r');
hold on; scatter(dirD_b, log(narrow_b./wide_b),'b');
hold on; plot(dirD, y);
hold on; plot(0:180, zeros(181,1),'k:');
legend('Monkey A','Monkey B','linear regression'); 
xlabel('delta theta'); ylabel('FR log ratio'); xlim([0 180]);
title(['Correlation between FR log ratio and delta-theta, rho =  ' num2str(rho) ', p = ' num2str(p) '']); 

%% PSTH - categorized by delta-theta

clear; clc;
load('A.mat');
wide_a = wide; narrow_a = narrow;
load('B.mat');
wide_b = wide; narrow_b = narrow;
load('delta_theta.mat');

% categorize
dCrit = 45; % direction (degree) for categorization criterion
large_a = find(dCrit<dirD_a);
small_a = find(dCrit>dirD_a);
large_b = find(dCrit<dirD_b);
small_b = find(dCrit>dirD_b);

% PSTH
figure;
subplot(1,2,1)
x = -100:0;
errorbar(x,nanmean(wide_a(large_a,200:end)),nanstd(wide_a(large_a,200:end))/sqrt(size(wide_a(large_a,200:end),1)));
hold on;
errorbar(x,nanmean(narrow_a(large_a,200:end)), nanstd(narrow_a(large_a,200:end))/sqrt(size(narrow_a(large_a,200:end),1)));
xlabel('Time (ms)'); ylabel('Firing rate'); title('A, Fixation task, Large delta-theta');
legend('','','','large','','','','small');

subplot(1,2,2)
errorbar(x,nanmean(wide_a(small_a,200:end)),nanstd(wide_a(small_a,200:end))/sqrt(size(wide_a(small_a,200:end),1)));
hold on;
errorbar(x,nanmean(narrow_a(small_a,200:end)), nanstd(narrow_a(small_a,200:end))/sqrt(size(narrow_a(small_a,200:end),1)));
xlabel('Time (ms)'); ylabel('Firing rate'); title('A, Fixation task, Small delta-theta');

figure;
subplot(1,2,1)
x = -100:0;
errorbar(x,nanmean(wide_b(large_b,200:end)),nanstd(wide_b(large_b,200:end))/sqrt(size(wide_b(large_b,200:end),1)));
hold on;
errorbar(x,nanmean(narrow_b(large_b,200:end)), nanstd(narrow_b(large_b,200:end))/sqrt(size(narrow_b(large_b,200:end),1)));
xlabel('Time (ms)'); ylabel('Firing rate'); title('B, Fixation task, Large delta-theta');

subplot(1,2,2)
errorbar(x,nanmean(wide_b(small_b,200:end)),nanstd(wide_b(small_b,200:end))/sqrt(size(wide_b(small_b,200:end),1)));
hold on;
errorbar(x,nanmean(narrow_b(small_b,200:end)), nanstd(narrow_b(small_b,200:end))/sqrt(size(narrow_b(small_b,200:end),1)));
xlabel('Time (ms)'); ylabel('Firing rate'); title('B, Fixation task, Small delta-theta');

% t-test
ttestResult_large_b = [];
[h,p,ci,stats] = ttest(mean(wide_b(large_b,200:end),2),mean(narrow_b(large_b,200:end),2));
ttestResult_large_b(1,1) = h;
ttestResult_large_b(1,2) = p;
ttestResult_large_b(1,3) = stats.tstat;
ttestResult_large_b(1,4) = nanmean(mean(wide_b(large_b,200:end),2));
ttestResult_large_b(1,5) = nanmean(mean(narrow_b(large_b,200:end),2));

ttestResult_small_b = [];
[h,p,ci,stats] = ttest(mean(wide_b(small_b,200:end),2),mean(narrow_b(small_b,200:end),2));
ttestResult_small_b(1,1) = h;
ttestResult_small_b(1,2) = p;
ttestResult_small_b(1,3) = stats.tstat;
ttestResult_small_b(1,4) = nanmean(mean(wide_b(small_b,200:end),2));
ttestResult_small_b(1,5) = nanmean(mean(narrow_b(small_b,200:end),2));

ttestResult_large_a = [];
[h,p,ci,stats] = ttest(mean(wide_a(large_a,200:end),2),mean(narrow_a(large_a,200:end),2));
ttestResult_large_a(1,1) = h;
ttestResult_large_a(1,2) = p;
ttestResult_large_a(1,3) = stats.tstat;
ttestResult_large_a(1,4) = nanmean(mean(wide_a(large_a,200:end),2));
ttestResult_large_a(1,5) = nanmean(mean(narrow_a(large_a,200:end),2));

ttestResult_small_a = [];
[h,p,ci,stats] = ttest(mean(wide_a(small_a,200:end),2),mean(narrow_a(small_a,200:end),2));
ttestResult_small_a(1,1) = h;
ttestResult_small_a(1,2) = p;
ttestResult_small_a(1,3) = stats.tstat;
ttestResult_small_a(1,4) = nanmean(mean(wide_a(small_a,200:end),2));
ttestResult_small_a(1,5) = nanmean(mean(narrow_a(small_a,200:end),2));