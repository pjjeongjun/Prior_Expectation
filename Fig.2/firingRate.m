% Mean firing rate in the time window from spike latency to 100ms
% bar graphs in Fig. 2, a-c

clear; clc; close all;

% load spike data
load('spikes.mat');

%% Monkey A

% Firing rate
b1c100 = mSpikes_A_b1c100; % wide-prior, high contrast
b2c100 = mSpikes_A_b2c100; % narrow-prior, high contrast
b1c008 = mSpikes_A_b1c008; % wide-prior, low contrast
b2c008 = mSpikes_A_b2c008; % narrow-prior, low contrast

% Spike latency (high contrast)
delete_idx = idx2_H_A; % index for 'spk lat > 600ms'
spkLat = spkLat_H_A; % spk latency

% Wide, High
frb1c100 = [];
for i = 1:size(b1c100,1)
    if ~isempty(find(i == delete_idx))
        frb1c100 = [frb1c100; nan];
    else
        frb1c100 = [frb1c100; mean(b1c100(i, spkLat(i):spkLat(i)+99))];
    end
end

% Narrow, High
frb2c100 = [];
for i = 1:size(b2c100,1)
    if ~isempty(find(i == delete_idx))
        frb2c100 = [frb2c100; nan];
    else
        frb2c100 = [frb2c100; mean(b2c100(i, spkLat(i):spkLat(i)+99))];
    end
end

% ttest
ttestResult100 = [];
for sample = 1:size(frb1c100,2)
    [h,p,ci,stats] = ttest(frb1c100(:,sample),frb2c100(:,sample));
    ttestResult100(sample,1) = h;
    ttestResult100(sample,2) = p;
    ttestResult100(sample,3) = stats.tstat;
    ttestResult100(sample,4) = nanmean(frb1c100(:,sample));
    ttestResult100(sample,5) = nanmean(frb2c100(:,sample));
end
save('ttest_high_A.mat','ttestResult100');

fig = figure;
c = categorical({'wide-prior','narrow-prior'});
mm = [nanmean(frb1c100),nanmean(frb2c100)];
sd = [nanstd(frb1c100)/sqrt(size(frb1c100,1)),nanstd(frb2c100)/sqrt(size(frb2c100,1))];
errorbar(c,mm,sd);
hold on;
bar(c(1),mm(1),'k');
hold on;
bar(c(2),mm(2),'r');
title(['High contrast, p = ' num2str(ttestResult100(1,2)) '']);
ylim([0 60]); ylabel('Firing rate');
savefig('bar_high_A.fig')

% Spike latency (low contrast)
delete_idx = idx2_L_A; % index for 'spk lat > 600ms'
spkLat = spkLat_L_A; % spk latency

% Wide, Low
frb1c008 = [];
for i = 1:size(b1c008,1)
    if ~isempty(find(i == delete_idx))
        frb1c008 = [frb1c008; nan];
    else
        frb1c008 = [frb1c008; mean(b1c008(i, spkLat(i):spkLat(i)+99))];
    end
end

% Narrow, Low
frb2c008 = [];
for i = 1:size(b2c008,1)
    if ~isempty(find(i == delete_idx))
        frb2c008 = [frb2c008; nan];
    else
        frb2c008 = [frb2c008; mean(b2c008(i, spkLat(i):spkLat(i)+99))];
    end
end

% ttest
ttestResult008 = [];
for sample = 1:size(frb1c008,2)
    [h,p,ci,stats] = ttest(frb1c008(:,sample),frb2c008(:,sample));
    ttestResult008(sample,1) = h;
    ttestResult008(sample,2) = p;
    ttestResult008(sample,3) = stats.tstat;
    ttestResult008(sample,4) = mean(frb1c008(:,sample));
    ttestResult008(sample,5) = mean(frb2c008(:,sample));
end
save('ttest_low_A.mat','ttestResult008');

fig = figure;
c = categorical({'wide-prior','narrow-prior'});
mm = [nanmean(frb1c008),nanmean(frb2c008)];
sd = [nanstd(frb1c008)/sqrt(size(frb1c008,1)),nanstd(frb2c008)/sqrt(size(frb2c008,1))];
errorbar(c,mm,sd);
hold on;
bar(c(1),mm(1),'k')
hold on;
bar(c(2),mm(2),'r')
title(['Low contrast, p = ' num2str(ttestResult008(1,2)) '']);
ylim([0 60]); ylabel('Firing rate');
savefig('bar_low_A.fig')

%% Monkey B

% Firing rate
b1c100 = mSpikes_B_b1c100; % wide-prior, high contrast
b2c100 = mSpikes_B_b2c100; % narrow-prior, high contrast
b1c008 = mSpikes_B_b1c008; % wide-prior, low contrast
b2c008 = mSpikes_B_b2c008; % narrow-prior, low contrast

% Spike latency (high contrast)
delete_idx = idx2_H_B; % index for 'spk lat > 600ms'
spkLat = spkLat_H_B; % spk latency

% Wide, High
frb1c100 = [];
for i = 1:size(b1c100,1)
    if ~isempty(find(i == delete_idx))
        frb1c100 = [frb1c100; nan];
    else
        frb1c100 = [frb1c100; mean(b1c100(i, spkLat(i):spkLat(i)+99))];
    end
end

% Narrow, High
frb2c100 = [];
for i = 1:size(b2c100,1)
    if ~isempty(find(i == delete_idx))
        frb2c100 = [frb2c100; nan];
    else
        frb2c100 = [frb2c100; mean(b2c100(i, spkLat(i):spkLat(i)+99))];
    end
end

% ttest
ttestResult100 = [];
for sample = 1:size(frb1c100,2)
    [h,p,ci,stats] = ttest(frb1c100(:,sample),frb2c100(:,sample));
    ttestResult100(sample,1) = h;
    ttestResult100(sample,2) = p;
    ttestResult100(sample,3) = stats.tstat;
    ttestResult100(sample,4) = mean(frb1c100(:,sample));
    ttestResult100(sample,5) = mean(frb2c100(:,sample));
end
save('ttest_high_B.mat','ttestResult100');

fig = figure;
c = categorical({'wide-prior','narrow-prior'});
mm = [nanmean(frb1c100),nanmean(frb2c100)];
sd = [nanstd(frb1c100)/sqrt(size(frb1c100,1)),nanstd(frb2c100)/sqrt(size(frb2c100,1))];
errorbar(c,mm,sd);
hold on;
bar(c(1),mm(1),'k')
hold on;
bar(c(2),mm(2),'r')
title(['High contrast, p = ' num2str(ttestResult100(1,2)) '']);
ylim([0 60]); ylabel('Firing rate');
savefig('bar_high_B.fig')

% Spike latency (low contrast)
delete_idx = idx2_L_B; % index for 'spk lat > 600ms'
spkLat = spkLat_L_B; % spk latency

% Wide, Low
frb1c008 = [];
for i = 1:size(b1c008,1)
    if ~isempty(find(i == delete_idx))
        frb1c008 = [frb1c008; nan];
    else
        frb1c008 = [frb1c008; mean(b1c008(i, spkLat(i):spkLat(i)+99))];
    end
end

% Narrow, Low
frb2c008 = [];
for i = 1:size(b2c008,1)
    if ~isempty(find(i == delete_idx))
        frb2c008 = [frb2c008; nan];
    else
        frb2c008 = [frb2c008; mean(b2c008(i, spkLat(i):spkLat(i)+99))];
    end
end

% t-test
ttestResult008 = [];
for sample = 1:size(frb1c008,2)
    [h,p,ci,stats] = ttest(frb1c008(:,sample),frb2c008(:,sample));
    ttestResult008(sample,1) = h;
    ttestResult008(sample,2) = p;
    ttestResult008(sample,3) = stats.tstat;
    ttestResult008(sample,4) = mean(frb1c008(:,sample));
    ttestResult008(sample,5) = mean(frb2c008(:,sample));
end
save('ttest_low_B.mat','ttestResult008');

fig = figure;
c = categorical({'wide-prior','narrow-prior'});
mm = [nanmean(frb1c008),nanmean(frb2c008)];
sd = [nanstd(frb1c008)/sqrt(size(frb1c008,1)),nanstd(frb2c008)/sqrt(size(frb2c008,1))];
errorbar(c,mm,sd);
hold on;
bar(c(1),mm(1),'k')
hold on;
bar(c(2),mm(2),'r')
title(['Low contrast, p = ' num2str(ttestResult008(1,2)) '']);
ylim([0 60]); ylabel('Firing rate');
savefig('bar_low_B.fig')

%% Monkey A, Low contrast, Sig/Non-sig behavioral effect days
% Firing rate
b1c100 = mSpikes_A_b1c100; % wide-prior, high contrast
b2c100 = mSpikes_A_b2c100; % narrow-prior, high contrast
b1c008 = mSpikes_A_b1c008; % wide-prior, low contrast
b2c008 = mSpikes_A_b2c008; % narrow-prior, low contrast

% Spike latency (low contrast)
delete_idx = idx2_L_A; % index for 'spk lat > 600ms'
spkLat = spkLat_L_A; % spk latency

% Index for 'Sig/Non-sig days' (behavioral effect)
highidxB = highidxB_A; % high contrast, not significant
highidxG = highidxG_A; % high contrast, significant
lowidxB = lowidxB_A; % low contrast, not significant
lowidxG = lowidxG_A; % low contrast, significant

day='Sig'; %'Sig' for sig days / 'Non' for non-sig days
if strcmp(day,'Sig')
    didx=lowidxG; 
elseif strcmp(day,'Non')
    didx=lowidxB;
end

% wide-prior, low contrast
frb1c008 = [];
for i = 1:length(didx)
    if ~isempty(find(didx(i) == delete_idx))
        frb1c008 = [frb1c008; nan];
    else
        frb1c008 = [frb1c008; mean(b1c008(didx(i), spkLat(didx(i)):spkLat(didx(i))+99))];
    end
end

% narrow-prior, low contrast
frb2c008 = [];
for i = 1:length(didx)
    if ~isempty(find(didx(i) == delete_idx))
        frb2c008 = [frb2c008; nan];
    else
        frb2c008 = [frb2c008; mean(b2c008(didx(i), spkLat(didx(i)):spkLat(didx(i))+99))];
    end
end

% t-test
ttestResult008 = [];
for sample = 1:size(frb1c008,2)
    [h,p,ci,stats] = ttest(frb1c008(:,sample),frb2c008(:,sample));
    ttestResult008(sample,1) = h;
    ttestResult008(sample,2) = p;
    ttestResult008(sample,3) = stats.tstat;
    ttestResult008(sample,4) = nanmean(frb1c008(:,sample));
    ttestResult008(sample,5) = nanmean(frb2c008(:,sample));
end
save(['ttest_low_' day '_A.mat'],'ttestResult008');

fig = figure;
c = categorical({'wide-prior','narrow-prior'});
mm = [nanmean(frb1c008),nanmean(frb2c008)];
sd = [nanstd(frb1c008)/sqrt(size(frb1c008,1)),nanstd(frb2c008)/sqrt(size(frb2c008,1))];
errorbar(c,mm,sd);
hold on;
bar(c(1),mm(1),'k')
hold on;
bar(c(2),mm(2),'r')
title(['' day ', Low contrast, p = ' num2str(ttestResult008(1,2)) '']);
ylim([0 60]); ylabel('Firing rate');
savefig(['bar_low_' day '_A.fig'])

save(['data_' day '_A.mat'],'frb1c008','frb2c008');

%% Monkey B, Low, Sig/Non-sig behavioral effect days
% Firing rate
b1c100 = mSpikes_B_b1c100; % wide-prior, high contrast
b2c100 = mSpikes_B_b2c100; % narrow-prior, high contrast
b1c008 = mSpikes_B_b1c008; % wide-prior, low contrast
b2c008 = mSpikes_B_b2c008; % narrow-prior, low contrast

% Spike latency (low contrast)
delete_idx = idx2_L_B; % index for 'spk lat > 600ms'
spkLat = spkLat_L_B; % spk latency

% Index for 'Sig/Non sig days' (behavioral effect)
highidxB = highidxB_B; % high contrast, not significant
highidxG = highidxG_B; % high contrast, significant
lowidxB = lowidxB_B; % low contrast, not significant
lowidxG = lowidxG_B; % low contrast, significant

day='Sig'; %'Sig' for sig days / 'Non' for non-sig days
if strcmp(day,'Sig')
    didx=lowidxG; 
elseif strcmp(day,'Non')
    didx=lowidxB;
end

% wide-prior, low contrast
frb1c008 = [];
for i = 1:length(didx)
    if ~isempty(find(didx(i) == delete_idx))
        frb1c008 = [frb1c008; nan];
    else
        frb1c008 = [frb1c008; mean(b1c008(didx(i), spkLat(didx(i)):spkLat(didx(i))+99))];
    end
end

% narrow-prior, low contrast
frb2c008 = [];
for i = 1:length(didx)
    if ~isempty(find(didx(i) == delete_idx))
        frb2c008 = [frb2c008; nan];
    else
        frb2c008 = [frb2c008; mean(b2c008(didx(i), spkLat(didx(i)):spkLat(didx(i))+99))];
    end
end

% t-test
ttestResult008 = [];
for sample = 1:size(frb1c008,2)
    [h,p,ci,stats] = ttest(frb1c008(:,sample),frb2c008(:,sample));
    ttestResult008(sample,1) = h;
    ttestResult008(sample,2) = p;
    ttestResult008(sample,3) = stats.tstat;
    ttestResult008(sample,4) = nanmean(frb1c008(:,sample));
    ttestResult008(sample,5) = nanmean(frb2c008(:,sample));
end
save(['ttest_low_' day '_B.mat'],'ttestResult008');

fig = figure;
c = categorical({'wide-prior','narrow-prior'});
mm = [nanmean(frb1c008),nanmean(frb2c008)];
sd = [nanstd(frb1c008)/sqrt(size(frb1c008,1)),nanstd(frb2c008)/sqrt(size(frb2c008,1))];
errorbar(c,mm,sd);
hold on;
bar(c(1),mm(1),'k')
hold on;
bar(c(2),mm(2),'r')
title(['' day ', Low contrast, p = ' num2str(ttestResult008(1,2)) '']);
ylim([0 60]); ylabel('Firing rate');
savefig(['bar_low_' day '_B.fig'])

save(['data_' day '_B.mat'],'frb1c008','frb2c008');

%% Combined, Low, Sig/Non-sig behavioral effect days
load('data_Sig_A.mat');
agfrb1c008=frb1c008;
agfrb2c008=frb2c008;
load('data_Non_A.mat');
abfrb1c008=frb1c008;
abfrb2c008=frb2c008;
load('data_Sig_B.mat');
bgfrb1c008=frb1c008;
bgfrb2c008=frb2c008;
load('data_Non_B.mat');
bbfrb1c008=frb1c008;
bbfrb2c008=frb2c008;

%sig
gfrb1c008=[agfrb1c008; bgfrb1c008];
gfrb2c008=[agfrb2c008; bgfrb2c008];
%non-sig
bfrb1c008=[abfrb1c008; bbfrb1c008];
bfrb2c008=[abfrb2c008; bbfrb2c008];

% ttest - sig days
ttestResult008 = [];
for sample = 1:size(gfrb1c008,2)
    [h,p,ci,stats] = ttest(gfrb1c008(:,sample),gfrb2c008(:,sample));
    ttestResult008(sample,1) = h;
    ttestResult008(sample,2) = p;
    ttestResult008(sample,3) = stats.tstat;
    ttestResult008(sample,4) = nanmean(gfrb1c008(:,sample));
    ttestResult008(sample,5) = nanmean(gfrb2c008(:,sample));
end
save('ttest_low_Sig.mat','ttestResult008');

fig = figure;
c = categorical({'wide-prior','narrow-prior'});
mm = [nanmean(gfrb1c008),nanmean(gfrb2c008)];
sd = [nanstd(gfrb1c008)/sqrt(size(gfrb1c008,1)),nanstd(gfrb2c008)/sqrt(size(gfrb2c008,1))];
errorbar(c,mm,sd);
hold on;
bar(c(1),mm(1),'k')
hold on;
bar(c(2),mm(2),'r')
title(['Low contrast, p = ' num2str(ttestResult008(1,2)) '']);
ylim([0 60]); ylabel('Firing rate');
savefig('bar_low_sig.fig')

% ttest - non-sig days
ttestResult008 = [];
for sample = 1:size(bfrb1c008,2)
    [h,p,ci,stats] = ttest(bfrb1c008(:,sample),bfrb2c008(:,sample));
    ttestResult008(sample,1) = h;
    ttestResult008(sample,2) = p;
    ttestResult008(sample,3) = stats.tstat;
    ttestResult008(sample,4) = nanmean(bfrb1c008(:,sample));
    ttestResult008(sample,5) = nanmean(bfrb2c008(:,sample));
end
save('ttest_low_B.mat','ttestResult008');

fig = figure;
c = categorical({'wide-prior','narrow-prior'});
mm = [nanmean(bfrb1c008),nanmean(bfrb2c008)];
sd = [nanstd(bfrb1c008)/sqrt(size(bfrb1c008,1)),nanstd(bfrb2c008)/sqrt(size(bfrb2c008,1))];
errorbar(c,mm,sd);
hold on;
bar(c(1),mm(1),'k')
hold on;
bar(c(2),mm(2),'r')
title(['Low contrast, p = ' num2str(ttestResult008(1,2)) '']);
ylim([0 60]); ylabel('Firing rate');
savefig('bar_low_B.fig')

%difference - sig vs non-sig days
diffbad=bbfrb1c008-bbfrb2c008;
diffgood=bgfrb1c008-bgfrb2c008;
ttestResult = [];
[h,p,ci,stats] = ttest2(diffbad(:,sample),diffgood(:,sample));
ttestResult(sample,1) = h;
ttestResult(sample,2) = p;
ttestResult(sample,3) = stats.tstat;
ttestResult(sample,4) = nanmean(diffbad(:,sample));
ttestResult(sample,5) = nanmean(diffgood(:,sample));