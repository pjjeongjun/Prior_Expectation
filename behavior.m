% direction, speed, latency of smooth pursuit eye movements for prior direction
% Fig. 1 d,e and fig. S1 c-f

clear; clc; close all;

% load data
load('behavior.mat');

%% direction
direction = direction_day_prior_A; 
% direction_prior_day_A: prior direction, monkey A, for each day
% direction_prior_day_B: prior direction, monkey B, for each day
% direction_prior_cell_A: prior direction, monkey A, for each cell
% direction_prior_cell_B: prior direction, monkey B, for each cell

trCrit = 50;

% mean
db1c100mean = []; db2c100mean = []; db1c008mean = []; db2c008mean = [];
for i = 1:size(direction,1)
    if length(direction{i,1}) >= trCrit && length(direction{i,2}) >= trCrit && length(direction{i,3}) >= trCrit && length(direction{i,4}) >= trCrit
    db1c100mean = [db1c100mean; mean(direction{i,1})]; % wide-prior, high contrast
    db2c100mean = [db2c100mean; mean(direction{i,2})]; % narrow-prior, high contrast
    db1c008mean = [db1c008mean; mean(direction{i,3})]; % wide-prior, low contrast
    db2c008mean = [db2c008mean; mean(direction{i,4})]; % narrow-prior, low contrast
    end
end

% figure: direction mean
fig = figure('position', [0, 40, 900, 450]);
c = categorical({'wide-prior','narrow-prior'});
subplot(1,2,1)
bar(c(1),mean(db1c100mean));
hold on;
errorbar(c(1),mean(db1c100mean),std(db1c100mean)/sqrt(size(direction,1)))
hold on;
bar(c(2),mean(db2c100mean));
hold on;
errorbar(c(2),mean(db2c100mean),std(db2c100mean)/sqrt(size(direction,1)))
title('high contrast');

subplot(1,2,2)
bar(c(1),mean(db1c008mean));
hold on;
errorbar(c(1),mean(db1c008mean),std(db1c008mean)/sqrt(size(direction,1)))
hold on;
bar(c(2),mean(db2c008mean));
hold on;
errorbar(c(2),mean(db2c008mean),std(db2c008mean)/sqrt(size(direction,1)))
title('low contrast');

savefig(['mean_direction_trCrit' num2str(trCrit) '.fig']);
saveas(fig, ['mean_direction_trCrit' num2str(trCrit) '.png']);

% ttest
ttestResult100 = [];
[h,p,ci,stats] = ttest(db1c100mean,db2c100mean);
ttestResult100(1,1) = h;
ttestResult100(1,2) = p;
ttestResult100(1,3) = stats.tstat;
ttestResult100(1,4) = nanmean(db1c100mean);
ttestResult100(1,5) = nanmean(db2c100mean);

ttestResult008 = [];
[h,p,ci,stats] = ttest(db1c008mean,db2c008mean);
ttestResult008(1,1) = h;
ttestResult008(1,2) = p;
ttestResult008(1,3) = stats.tstat;
ttestResult008(1,4) = nanmean(db1c008mean);
ttestResult008(1,5) = nanmean(db2c008mean);

save(['ttest(mean_dir)_trCrit' num2str(trCrit) '.mat'],'ttestResult100','ttestResult008');

% std
db1c100std = []; db2c100std = []; db1c008std = []; db2c008std = [];
for i = 1:size(direction,1)
    if length(direction{i,1}) >= trCrit && length(direction{i,2}) >= trCrit && length(direction{i,3}) >= trCrit && length(direction{i,4}) >= trCrit
        db1c100std = [db1c100std; std(direction{i,1})]; % wide-prior, high contrast
        db2c100std = [db2c100std; std(direction{i,2})]; % narrow-prior, high contrast
        db1c008std = [db1c008std; std(direction{i,3})]; % wide-prior, low contrast
        db2c008std = [db2c008std; std(direction{i,4})]; % narrow-prior, low contrast
    end
end

save(['dirs_trCrit' num2str(trCrit) '.mat'],'db1c100mean','db2c100mean','db1c008mean','db2c008mean','db1c100std','db2c100std','db1c008std','db2c008std');

% figure: direction std
fig = figure('position', [0, 40, 900, 450]);
c = categorical({'wide-prior','narrow-prior'});
subplot(1,2,1)
bar(c(1),mean(db1c100std));
hold on;
errorbar(c(1),mean(db1c100std),std(db1c100std)/sqrt(size(direction,1)))
hold on;
bar(c(2),mean(db2c100std));
hold on;
errorbar(c(2),mean(db2c100std),std(db2c100std)/sqrt(size(direction,1)))
title('high contrast');

subplot(1,2,2)
bar(c(1),mean(db1c008std));
hold on;
errorbar(c(1),mean(db1c008std),std(db1c008std)/sqrt(size(direction,1)))
hold on;
bar(c(2),mean(db2c008std));
hold on;
errorbar(c(2),mean(db2c008std),std(db2c008std)/sqrt(size(direction,1)))
title('low contrast');

savefig(['std_direction_trCrit' num2str(trCrit) '.fig']);

%ttest
ttestResult100 = [];
[h,p,ci,stats] = ttest(db1c100std,db2c100std);
ttestResult100(1,1) = h;
ttestResult100(1,2) = p;
ttestResult100(1,3) = stats.tstat;
ttestResult100(1,4) = nanmean(db1c100std);
ttestResult100(1,5) = nanmean(db2c100std);

ttestResult008 = [];
[h,p,ci,stats] = ttest(db1c008std,db2c008std);
ttestResult008(1,1) = h;
ttestResult008(1,2) = p;
ttestResult008(1,3) = stats.tstat;
ttestResult008(1,4) = nanmean(db1c008std);
ttestResult008(1,5) = nanmean(db2c008std);

save(['ttest(std_dir)_trCrit' num2str(trCrit) '.mat'],'ttestResult100','ttestResult008');
save('std.mat','db1c100std','db2c100std','db1c008std','db2c008std');

% difference between high and low
load('std.mat');
high = db1c100std-db2c100std;
low = db1c008std-db2c008std;

ttestResult = [];
[h,p,ci,stats] = ttest(high,low);
ttestResult(1,1) = h;
ttestResult(1,2) = p;
ttestResult(1,3) = stats.tstat;
ttestResult(1,4) = nanmean(high);
ttestResult(1,5) = nanmean(low);

save('ttest_contrast.mat','ttestResult');

%% speed
speed = speed_day_prior_A; 
%speed_day_prior_A: prior direction, monkey A, for each day
%speed_day_prior_B: prior direction, monkey B, for each day
%speed_cell_prior_A: prior direction, monkey A, for each cell
%speed_cell_prior_B: prior direction, monkey B, for each cell

trCrit = 50;

% mean
sb1c100mean = []; sb2c100mean = []; sb1c008mean = []; sb2c008mean = [];
for i = 1:size(speed,1)
    if length(speed{i,1}) >= trCrit && length(speed{i,2}) >= trCrit && length(speed{i,3}) >= trCrit && length(speed{i,4}) >= trCrit
        sb1c100mean = [sb1c100mean; mean(speed{i,1})]; % wide-prior, high contrast
        sb2c100mean = [sb2c100mean; mean(speed{i,2})]; % narrow-prior, high contrast
        sb1c008mean = [sb1c008mean; mean(speed{i,3})]; % wide-prior, low contrast
        sb2c008mean = [sb2c008mean; mean(speed{i,4})]; % narrow-prior, low contrast
    end
end

% figure: speed mean
fig = figure('position', [0, 40, 900, 450]);
c = categorical({'wide-prior','narrow-prior'});
subplot(1,2,1)
bar(c(1),mean(sb1c100mean));
hold on;
errorbar(c(1),mean(sb1c100mean),std(sb1c100mean)/size(speed,1))
hold on;
bar(c(2),mean(sb2c100mean));
hold on;
errorbar(c(2),mean(sb2c100mean),std(sb2c100mean)/size(speed,1))
title('high contrast');

subplot(1,2,2)
bar(c(1),mean(sb1c008mean));
hold on;
errorbar(c(1),mean(sb1c008mean),std(sb1c008mean)/size(speed,1))
hold on;
bar(c(2),mean(sb2c008mean));
hold on;
errorbar(c(2),mean(sb2c008mean),std(sb2c008mean)/size(speed,1))
title('low contrast');

savefig(['mean_speed_trCrit' num2str(trCrit) '.fig']);

%ttest
ttestResult100 = [];
[h,p,ci,stats] = ttest(sb1c100mean,sb2c100mean);
ttestResult100(1,1) = h;
ttestResult100(1,2) = p;
ttestResult100(1,3) = stats.tstat;
ttestResult100(1,4) = nanmean(sb1c100mean);
ttestResult100(1,5) = nanmean(sb2c100mean);

ttestResult008 = [];
[h,p,ci,stats] = ttest(sb1c008mean,sb2c008mean);
ttestResult008(1,1) = h;
ttestResult008(1,2) = p;
ttestResult008(1,3) = stats.tstat;
ttestResult008(1,4) = nanmean(sb1c008mean);
ttestResult008(1,5) = nanmean(sb2c008mean);

save(['ttest(mean_spd)_trCrit' num2str(trCrit) '.mat'],'ttestResult100','ttestResult008');

% std
sb1c100std = []; sb2c100std = []; sb1c008std = []; sb2c008std = [];
for i = 1:size(speed,1)
    if length(speed{i,1}) >= trCrit && length(speed{i,2}) >= trCrit && length(speed{i,3}) >= trCrit && length(speed{i,4}) >= trCrit
        sb1c100std = [sb1c100std; std(speed{i,1})]; % wide-prior, high contrast
        sb2c100std = [sb2c100std; std(speed{i,2})]; % narrow-prior, high contrast
        sb1c008std = [sb1c008std; std(speed{i,3})]; % wide-prior, low contrast
        sb2c008std = [sb2c008std; std(speed{i,4})]; % narrow-prior, low contrast
    end
end

save(['spds_trCrit' num2str(trCrit) '.mat'], 'sb1c100mean','sb2c100mean','sb1c008mean','sb2c008mean','sb1c100std','sb2c100std','sb1c008std','sb2c008std');

% figure: speed std
fig = figure('position', [0, 40, 900, 450]);
c = categorical({'wide-prior','narrow-prior'});
subplot(1,2,1)
bar(c(1),mean(sb1c100std));
hold on;
errorbar(c(1),mean(sb1c100std),std(sb1c100std)/size(speed,1))
hold on;
bar(c(2),mean(sb2c100std));
hold on;
errorbar(c(2),mean(sb2c100std),std(sb2c100std)/size(speed,1))
title('high contrast');

subplot(1,2,2)
bar(c(1),mean(sb1c008std));
hold on;
errorbar(c(1),mean(sb1c008std),std(sb1c008std)/size(speed,1))
hold on;
bar(c(2),mean(sb2c008std));
hold on;
errorbar(c(2),mean(sb2c008std),std(sb2c008std)/size(speed,1))
title('low contrast');

savefig(['std_speed_trCrit' num2str(trCrit) '.fig']);

% ttest
ttestResult100 = [];
[h,p,ci,stats] = ttest(sb1c100std,sb2c100std);
ttestResult100(1,1) = h;
ttestResult100(1,2) = p;
ttestResult100(1,3) = stats.tstat;
ttestResult100(1,4) = nanmean(sb1c100std);
ttestResult100(1,5) = nanmean(sb2c100std);

ttestResult008 = [];
[h,p,ci,stats] = ttest(sb1c008std,sb2c008std);
ttestResult008(1,1) = h;
ttestResult008(1,2) = p;
ttestResult008(1,3) = stats.tstat;
ttestResult008(1,4) = nanmean(sb1c008std);
ttestResult008(1,5) = nanmean(sb2c008std);

save(['ttest(std_spd)_trCrit' num2str(trCrit) '.mat'],'ttestResult100','ttestResult008');

%% latency
latency = latency_day_prior_A; 
%latency_day_prior_A: prior direction, monkey A, for each day
%latency_day_prior_B: prior direction, monkey B, for each day
%latency_cell_prior_A: prior direction, monkey A, for each cell
%latency_cell_prior_B: prior direction, monkey B, for each cell

trCrit = 50;

% mean
lb1c100mean = []; lb2c100mean = []; lb1c008mean = []; lb2c008mean = [];
for i = 1:size(latency,1)
    if length(latency{i,1}) >= trCrit && length(latency{i,2}) >= trCrit && length(latency{i,3}) >= trCrit && length(latency{i,4}) >= trCrit
        lb1c100mean = [lb1c100mean; mean(latency{i,1})]; % wide-prior, high contrast
        lb2c100mean = [lb2c100mean; mean(latency{i,2})]; % narrow-prior, high contrast
        lb1c008mean = [lb1c008mean; mean(latency{i,3})]; % wide-prior, low contrast
        lb2c008mean = [lb2c008mean; mean(latency{i,4})]; % narrow-prior, low contrast
    end
end

% figure: latency mean
fig = figure('position', [0, 40, 900, 450]);
c = categorical({'block1','block2'});
subplot(1,2,1)
bar(c(1),mean(lb1c100mean));
hold on;
errorbar(c(1),mean(lb1c100mean),std(lb1c100mean)/size(latency,1))
hold on;
bar(c(2),mean(lb2c100mean));
hold on;
errorbar(c(2),mean(lb2c100mean),std(lb2c100mean)/size(latency,1))
title('high contrast');

subplot(1,2,2)
bar(c(1),mean(lb1c008mean));
hold on;
errorbar(c(1),mean(lb1c008mean),std(lb1c008mean)/size(latency,1))
hold on;
bar(c(2),mean(lb2c008mean));
hold on;
errorbar(c(2),mean(lb2c008mean),std(lb2c008mean)/size(latency,1))
title('low contrast');

savefig(['mean_latency_trCrit' num2str(trCrit) '.fig']);

% ttest
ttestResult100 = [];
[h,p,ci,stats] = ttest(lb1c100mean,lb2c100mean);
ttestResult100(1,1) = h;
ttestResult100(1,2) = p;
ttestResult100(1,3) = stats.tstat;
ttestResult100(1,4) = nanmean(lb1c100mean);
ttestResult100(1,5) = nanmean(lb2c100mean);

ttestResult008 = [];
[h,p,ci,stats] = ttest(lb1c008mean,lb2c008mean);
ttestResult008(1,1) = h;
ttestResult008(1,2) = p;
ttestResult008(1,3) = stats.tstat;
ttestResult008(1,4) = nanmean(lb1c008mean);
ttestResult008(1,5) = nanmean(lb2c008mean);

save(['ttest(mean_lat)_trCrit' num2str(trCrit) '.mat'],'ttestResult100','ttestResult008');

% std
lb1c100std = []; lb2c100std = []; lb1c008std = []; lb2c008std = [];
for i = 1:size(latency,1)
    if length(latency{i,1}) >= trCrit && length(latency{i,2}) >= trCrit && length(latency{i,3}) >= trCrit && length(latency{i,4}) >= trCrit
        lb1c100std = [lb1c100std; std(latency{i,1})]; % wide-prior, high contrast
        lb2c100std = [lb2c100std; std(latency{i,2})]; % narrow-prior, high contrast
        lb1c008std = [lb1c008std; std(latency{i,3})]; % wide-prior, low contrast
        lb2c008std = [lb2c008std; std(latency{i,4})]; % narrow-prior, low contrast
    end
end

save(['lats_trCrit' num2str(trCrit) '.mat'],'lb1c100mean','lb2c100mean','lb1c008mean','lb2c008mean','lb1c100std','lb2c100std','lb1c008std','lb2c008std');

% figure: latency std
fig = figure('position', [0, 40, 900, 450]);
c = categorical({'block1','block2'});
subplot(1,2,1)
bar(c(1),mean(lb1c100std));
hold on;
errorbar(c(1),mean(lb1c100std),std(lb1c100std)/size(latency,1))
hold on;
bar(c(2),mean(lb2c100std));
hold on;
errorbar(c(2),mean(lb2c100std),std(lb2c100std)/size(latency,1))
title('high contrast');

subplot(1,2,2)
bar(c(1),mean(lb1c008std));
hold on;
errorbar(c(1),mean(lb1c008std),std(lb1c008std)/size(latency,1))
hold on;
bar(c(2),mean(lb2c008std));
hold on;
errorbar(c(2),mean(lb2c008std),std(lb2c008std)/size(latency,1))
title('low contrast');

savefig(['std_latency_trCrit' num2str(trCrit) '.fig']);

% ttest
ttestResult100 = [];
[h,p,ci,stats] = ttest(lb1c100std,lb2c100std);
ttestResult100(1,1) = h;
ttestResult100(1,2) = p;
ttestResult100(1,3) = stats.tstat;
ttestResult100(1,4) = nanmean(lb1c100std);
ttestResult100(1,5) = nanmean(lb2c100std);

ttestResult008 = [];
[h,p,ci,stats] = ttest(lb1c008std,lb2c008std);
ttestResult008(1,1) = h;
ttestResult008(1,2) = p;
ttestResult008(1,3) = stats.tstat;
ttestResult008(1,4) = nanmean(lb1c008std);
ttestResult008(1,5) = nanmean(lb2c008std);

save(['ttest(std_lat)_trCrit' num2str(trCrit) '.mat'],'ttestResult100','ttestResult008');