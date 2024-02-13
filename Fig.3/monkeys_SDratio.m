%% std ratio of pursuit directions
% Fig. 3, e-f

% load behavioral data
load('behavior.mat');

% behavioral result - Monkey A
direction = direction_day_prior_A;

ab1c100=[];ab2c100=[];ab1c008=[];ab2c008=[];
for a = 1:size(direction,1)
    ab1c100 = [ab1c100; std(direction{a,1})];
    ab2c100 = [ab2c100; std(direction{a,2})];
    ab1c008 = [ab1c008; std(direction{a,3})];
    ab2c008 = [ab2c008; std(direction{a,4})];
end

% behavioral result - Monkey B
direction = direction_day_prior_B;

bb1c100=[];bb2c100=[];bb1c008=[];bb2c008=[];
for b = 1:size(direction,1)
    bb1c100 = [bb1c100; std(direction{b,1})];
    bb2c100 = [bb2c100; std(direction{b,2})];
    bb1c008 = [bb1c008; std(direction{b,3})];
    bb2c008 = [bb2c008; std(direction{b,4})];
end

b1c100 = [ab1c100; bb1c100];
b2c100 = [ab2c100; bb2c100];
b1c008 = [ab1c008; bb1c008];
b2c008 = [ab2c008; bb2c008];

%std log ratio or ratio
% high = log(b2c100./b1c100);
% low = log(b2c008./b1c008);
high = b2c100./b1c100;
low = b2c008./b1c008;

%ttestH
ttestResultH = [];
[h,p,ci,stats] = ttest(high);
ttestResultH(1,1) = h;
ttestResultH(1,2) = p;
ttestResultH(1,3) = stats.tstat;
ttestResultH(1,4) = mean(high);

%ttestL
ttestResultL = [];
[h,p,ci,stats] = ttest(low);
ttestResultL(1,1) = h;
ttestResultL(1,2) = p;
ttestResultL(1,3) = stats.tstat;
ttestResultL(1,4) = mean(low);

save('stdratio.mat','high','low','ttestResultH','ttestResultL');
