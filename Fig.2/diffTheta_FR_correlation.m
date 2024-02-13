%% Correlation between FR ratio and difference of the directions(preferrred-target)
% Fig. 2, d-i

clear; clc; close all;

%load spike data
load('spikes.mat');

%load preferred and target directions
load('pf_tar_dirs.mat');

% Responses to prior direction 
% Monkey A
targetHighWide_A = mSpikes_A_b1c100; % wide-prior, high contrast
targetHighNarrow_A = mSpikes_A_b2c100; % narrow-prior, high contrast
targetLowWide_A = mSpikes_A_b1c008; % wide-prior, low contrast
targetLowNarrow_A = mSpikes_A_b2c008; % narrow-prior, low contrast
% Monkey B
targetHighWide_B = mSpikes_B_b1c100; % wide-prior, high contrast
targetHighNarrow_B = mSpikes_B_b2c100; % narrow-prior, high contrast
targetLowWide_B = mSpikes_B_b1c008; % wide-prior, low contrast
targetLowNarrow_B = mSpikes_B_b2c008; % narrow-prior, low contrast

% Preferred and Target directions
% Monkey A
pfDir_a = transpose(cell2mat(pfDir_A));
tarDir_a = transpose(cell2mat(tarDir_A));
% Monkey B
pfDir_b = transpose(cell2mat(pfDir_B));
tarDir_b = transpose(cell2mat(tarDir_B));

prolog = 400;
winDur = 80;
preW = 0;

dirD_a = abs(pfDir_a-tarDir_a);
tt1 = find(dirD_a > 180);
dirD_a(tt1) = 360-dirD_a(tt1);
dirD_b = abs(pfDir_b-tarDir_b);
tt2 = find(dirD_b > 180);
dirD_b(tt2) = 360-dirD_b(tt2);

latCrit = 500;
dirCrit = 180;

nIdxA_L = find(isnan(spkLat_L_A)| spkLat_L_A > latCrit | dirD_a > dirCrit);

tgW = targetLowWide_A;
tgN = targetLowNarrow_A;
latA = spkLat_L_A;
latA(nIdxA_L) = [];

tgW(nIdxA_L,:) = [];
tgN(nIdxA_L,:) = [];

nIdxB_L = find(isnan(spkLat_L_B) | spkLat_L_B > latCrit | dirD_b > dirCrit);

tgW = targetLowWide_B;
tgN = targetLowNarrow_B;
latB = spkLat_L_B;
latB(nIdxB_L) = [];

tgW(nIdxB_L,:) = [];
tgN(nIdxB_L,:) = [];

lowD_a = dirD_a;
lowD_b = dirD_b;
lowD_a(nIdxA_L) = [];
lowD_b(nIdxB_L) = [];

lowidx1 = find(~isnan(lowD_a));
lowidx2 = find(~isnan(lowD_b));

nIdxA_H = [];
nIdxA_H = find(isnan(spkLat_H_A)| spkLat_H_A > latCrit  | dirD_a > dirCrit);

tgW = targetHighWide_A;
tgN = targetHighNarrow_A;
latA = spkLat_H_A;
latA(nIdxA_H) = [];

tgW(nIdxA_H,:) = [];
tgN(nIdxA_H,:) = [];

nIdxB_H = [];
nIdxB_H = find(isnan(spkLat_H_B) | spkLat_H_B > latCrit | dirD_b > dirCrit);

tgW = targetHighWide_B;
tgN = targetHighNarrow_B;
latB = spkLat_H_B;
latB(nIdxB_H) = [];

tgW(nIdxB_H,:) = [];
tgN(nIdxB_H,:) = [];

highD_a = dirD_a;
highD_b = dirD_b;

highD_a(nIdxA_H) = [];
highD_b(nIdxB_H) = [];

highidx1 = find(~isnan(highD_a));
highidx2 = find(~isnan(highD_b));

stpWin = 60;

wStep = [0:5:200];

corrL = nan(length(wStep),2);
corrH = nan(length(wStep),2);

lowScatter = cell(length(wStep),1);
highScatter = cell(length(wStep),1);

corrL_A = nan(length(wStep),2);
corrH_A = nan(length(wStep),2);

lowScatter_A = cell(length(wStep),1);
highScatter_A = cell(length(wStep),1);

corrL_B = nan(length(wStep),2);
corrH_B = nan(length(wStep),2);

lowScatter_B = cell(length(wStep),1);
highScatter_B = cell(length(wStep),1);

for k = 1:length(wStep)
	rwin = [prolog-stpWin/2+wStep(k)+1:prolog+stpWin/2+wStep(k)];

	tgW = targetLowWide_A;
	tgN = targetLowNarrow_A;
	tgW(nIdxA_L,:) = [];
	tgN(nIdxA_L,:) = [];

	lowA_d = mean(tgN(:,rwin), 2) - mean(tgW(:, rwin), 2);
	lowA_r = mean(tgN(:, rwin), 2)./mean(tgW(:, rwin), 2);


	tgW = targetLowWide_B;
	tgN = targetLowNarrow_B;
	tgW(nIdxB_L,:) = [];
	tgN(nIdxB_L,:) = [];

	lowB_d = mean(tgN(:,rwin), 2) - mean(tgW(:, rwin), 2);
	lowB_r = mean(tgN(:, rwin), 2)./mean(tgW(:, rwin), 2);

	tgW = targetHighWide_A;
	tgN = targetHighNarrow_A;
	tgW(nIdxA_H,:) = [];
	tgN(nIdxA_H,:) = [];

	highA_d = mean(tgN(:,rwin), 2) - mean(tgW(:, rwin), 2);
	highA_r = mean(tgN(:, rwin), 2)./mean(tgW(:, rwin), 2);

	tgW = targetHighWide_B;
	tgN = targetHighNarrow_B;
	tgW(nIdxB_H,:) = [];
	tgN(nIdxB_H,:) = [];

	highB_d = mean(tgN(:,rwin), 2) - mean(tgW(:, rwin), 2);
	highB_r = mean(tgN(:, rwin), 2)./mean(tgW(:, rwin), 2);
	
	la_r = lowA_r(lowidx1);
	lb_r = lowB_r(lowidx2);
	ha_r = highA_r(highidx1);
	hb_r = highB_r(highidx2);
	
	la_d = lowD_a(lowidx1);
	lb_d = lowD_b(lowidx2);
	ha_d = highD_a(highidx1);
	hb_d = highD_b(highidx2);
	
	nal = isnan(la_r);
	nbl = isnan(lb_r);
	nah = isnan(ha_r);
	nbh = isnan(hb_r);
	
	la_r(nal) = [];
	lb_r(nbl) = [];
	ha_r(nah) = [];
	hb_r(nbh) = [];
	
	la_d(nal) = [];
	lb_d(nbl) = [];
	ha_d(nah) = [];
	hb_d(nbh) = [];
	
	lowComb = [[la_d;lb_d] [la_r;lb_r]];
	highComb = [[ha_d;hb_d] [ha_r;hb_r]];

	lowComb_a = [[la_d] [la_r]];
	highComb_a = [[ha_d] [ha_r]];

	lowComb_b = [[lb_d] [lb_r]];
	highComb_b = [[hb_d] [hb_r]];
	
	lowScatter{k} = lowComb;
	highScatter{k} = highComb;

	lowScatter_A{k} = lowComb_a;
	highScatter_A{k} = highComb_a;

	lowScatter_B{k} = lowComb_b;
	highScatter_B{k} = highComb_b;
	
	[corrL(k, 1), corrL(k, 2)] = corr(lowComb(:,1), log(lowComb(:,2)), 'type', 'spearman');
	[corrH(k, 1), corrH(k, 2)] = corr(highComb(:,1), log(highComb(:,2)), 'type', 'spearman');

	[corrL_A(k, 1), corrL_A(k, 2)] = corr(lowComb_a(:,1), log(lowComb_a(:,2)), 'type', 'spearman');
	[corrH_A(k, 1), corrH_A(k, 2)] = corr(highComb_a(:,1), log(highComb_a(:,2)), 'type', 'spearman');

	[corrL_B(k, 1), corrL_B(k, 2)] = corr(lowComb_b(:,1), log(lowComb_b(:,2)), 'type', 'spearman');
	[corrH_B(k, 1), corrH_B(k, 2)] = corr(highComb_b(:,1), log(highComb_b(:,2)), 'type', 'spearman');
end

save('correlation.mat','corrH','corrL','corrH_A','corrL_A','corrH_B','corrL_B','lowComb','highComb','lowComb_a','highComb_a','lowComb_b','highComb_b')
save('ratio.mat','highScatter','lowScatter','highScatter_A','lowScatter_A','highScatter_B','lowScatter_B');

%% figure
clear; clc;
wStep = [0:5:200];

%correlation
load('correlation.mat');
figure;
plot(wStep,corrH(:,1),'k') % Combined
hold on;
plot(wStep,corrH_A(:,1),'r') % Monkey A
hold on;
plot(wStep,corrH_B(:,1),'b') % Monkey B
legend('Combined','Monkey A','Monkey B');
xlabel('Time (ms)'); ylim([-0.4 0.4]); ylabel('rho');
title('correlation between firing rate ratio and delta-theta, High contrast');
savefig('corrH.fig')

figure;
plot(wStep,corrL(:,1),'k') % Combined
hold on;
plot(wStep,corrL_A(:,1),'r') % Monkey A
hold on;
plot(wStep,corrL_B(:,1),'b') % Monkey B
legend('Combined','Monkey A','Monkey B');
xlabel('Time (ms)'); ylim([-0.4 0.4]); ylabel('rho');
title('correlation between firing rate ratio and delta-theta, Low contrast');
savefig('corrL.fig')

%scatter(85 ms: 55-115 ms)
load('ratio.mat');
figure;
scatter(highScatter_A{18}(:,1),highScatter_A{18}(:,2),'r'); hold on;
scatter(highScatter_B{18}(:,1),highScatter_B{18}(:,2),'b');
xlabel('Time (ms)'); ylabel('Firing rate ratio');
title('55-115 ms, Firing rate ratio, High contrast');
savefig('scatterH.fig')

figure;
scatter(lowScatter_A{18}(:,1),lowScatter_A{18}(:,2),'r'); hold on;
scatter(lowScatter_B{18}(:,1),lowScatter_B{18}(:,2),'b'); hold on;
xlabel('Time (ms)'); ylabel('Firing rate ratio');
title('55-115 ms, Firing rate ratio, Low contrast');
savefig('scatterL.fig')

%regression
degree = [ones(length(highScatter{18}(:,1)),1) highScatter{18}(:,1)];
[bh,~,~,~,~] = regress(highScatter{18}(:,2),degree);
x = 0:180;
yh = bh(2)*x+bh(1);

degree = [ones(length(lowScatter{18}(:,1)),1) lowScatter{18}(:,1)];
[bl,~,~,~,~] = regress(lowScatter{18}(:,2),degree);
x = 0:180;
yl = bl(2)*x+bl(1);

%% regression
clear; clc; close all;
load('correlation.mat');
load('ratio.mat');

wStep = 0:5:200;
for m = 1:3
    rvH = []; rvL = []; ciH = []; ciL = [];
    for k = 1:length(wStep) %k=18: 85 ms (55-115 ms)
        if m == 1
            regH = highScatter{k};   %Combined
        elseif m == 2
            regH = highScatter_A{k}; % Monkey A
        else
            regH = highScatter_B{k}; % Monkey B
        end
        regH(:,2) = log(regH(:,2));
        regH(find(isinf(regH(:,2))),:) = [];

        if m == 1
            regL = lowScatter{k};   % Combined
        elseif m == 2
            regL = lowScatter_A{k}; % Monkey A
        else
            regL = lowScatter_B{k}; % Monkey B
        end
        regL(:,2) = log(regL(:,2));
        regL(find(isinf(regL(:,2))),:) = [];

        [h,hint] = regress(regH(:,2),[ones(size(regH,1),1) regH(:,1)]);
        [l,lint] = regress(regL(:,2),[ones(size(regL,1),1) regL(:,1)]);
        rvH = [rvH; h(2)]; ciH = [ciH; hint(2,:)];
        rvL = [rvL; l(2)]; ciL = [ciL; lint(2,:)];
    end

    if m == 1
        save('regression_combined.mat','rvH','rvL','ciH','ciL')
    elseif m == 2
        save('regression_A.mat','rvH','rvL','ciH','ciL')
    else
        save('regression_B.mat','rvH','rvL','ciH','ciL')
    end

    %figure
    %high
    figure;
    subplot(1,2,1);
    y1 = ciH(:,2)';
    y2 = ciH(:,1)';
    xconf = [wStep wStep(end:-1:1)];
    yconf = [y1 y2(end:-1:1)];

    p = fill(xconf,yconf,'blue');
    p.FaceColor = [0.8 0.8 1];
    p.EdgeColor = 'none';
    hold on; plot(wStep,rvH,'b');
    hold on; plot(wStep,zeros(length(wStep),1),'k:');
    ylim([-0.005 0.002]);
    title('High contrast'); xlabel('Time (ms)'); ylabel('Regression coefficient');

    %low
    subplot(1,2,2)
    y1 = ciL(:,2)';
    y2 = ciL(:,1)';
    xconf = [wStep wStep(end:-1:1)];
    yconf = [y1 y2(end:-1:1)];

    p = fill(xconf,yconf,'red');
    p.FaceColor = [1 0.8 0.8];
    p.EdgeColor = 'none';
    hold on; plot(wStep,rvL,'r');
    hold on; plot(wStep,zeros(length(wStep),1),'k:');
    ylim([-0.005 0.002]);
    title('Low contrast'); xlabel('Time (ms)'); ylabel('Regression coefficient');
    if m == 1
        savefig('regression_Combined.fig')
    elseif m == 2
        savefig('regression_A.fig')
    else
        savefig('regression_B.fig')
    end
end

%% FDR (false discovery rate)
clear; clc;

time = 0:5:200;

load('correlation.mat');
range=1:length(time)-4;
[~,~,~,qH] = fdr(corrH(range,2));
[~,~,~,qL] = fdr(corrL(range,2));

high=time(find(qH<0.05));
low=time(find(qL<0.05));