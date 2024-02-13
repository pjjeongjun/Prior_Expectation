%% Support Vector Machine
% Fig. 3d

clear; clc; close all;

% load experimental data
load('Monkeys_cir.mat');
Monkeys = Data;

% simulate the neural data

% gamma parameters of prefer and null response distribution
% high, prefer
hp = Monkeys(1).Tuning_pfnull(1,:);
nanidx = find(isnan(hp));
hp(nanidx) = [];
params1 = gamfit(hp); % shape and scale parameter
% high, null
hn = Monkeys(1).Tuning_pfnull(2,:);
nanidx = find(isnan(hn));
hn(nanidx) = [];
params2 = gamfit(hn); % shape and scale parameter
% low, prefer
lp = Monkeys(2).Tuning_pfnull(1,:);
nanidx = find(isnan(lp));
lp(nanidx) = [];
params3 = gamfit(lp); % shape and scale parameter
% low, null
ln = Monkeys(2).Tuning_pfnull(2,:);
nanidx = find(isnan(ln));
ln(nanidx) = [];
params4 = gamfit(ln); % shape and scale parameter

% gamma parameters of sigma
% high
hs = Monkeys(1).Sigma.*180.0/pi;
nanidx = find(isnan(hs));
hs(nanidx) = [];
params5 = gamfit(hs); % shape and scale parameter
% low
ls = Monkeys(2).Sigma.*180.0/pi;
nanidx = find(isnan(ls));
ls(nanidx) = [];
params6 = gamfit(ls); % shape and scale parameter

% gamma parameters of cparams
%high
hc1 = Monkeys(1).cparams(:,1); %param(1)
hc2 = Monkeys(1).cparams(:,2); %param(2)
hc3 = Monkeys(1).cparams(:,3); %param(3)
didx = unique([find(isnan(hc1)) find(isnan(hc2)) find(isnan(hc3))]);
hc1(didx) = []; hc2(didx) = []; hc3(didx) = [];
params7 = [gamfit(hc1); gamfit(hc2); gamfit(hc3)]; % shape and scale parameter for param(1)/(2)/(3)
%low
lc1 = Monkeys(2).cparams(:,1); %param(1)
lc2 = Monkeys(2).cparams(:,2); %param(2)
lc3 = Monkeys(2).cparams(:,3); %param(3)
didx = unique([find(isnan(lc1)) find(isnan(lc2)) find(isnan(lc3))]);
lc1(didx) = []; lc2(didx) = []; lc3(didx) = [];
params8 = [gamfit(lc1); gamfit(lc2); gamfit(lc3)]; % shape and scale parameter for param(1)/(2)/(3)

save('gamparams.mat','params1','params2','params3','params4','params5','params6','params7','params8','hc1','hc2','hc3','lc1','lc2','lc3');

% make simulated neurons
for c = 1:2
    
    if c == 1
        contrast = 'h';
    elseif c == 2
        contrast = 'l';
    end
    
    %------set the parameters ------%
    twin = 18; %18 = 85 ms (55-115 ms)
    nMulti = 10; % number of neurons for each direction
    Npd = 360; % number of prefer directions
    numTrials = 200; % number of trials (make it similar to acutual trial num)
    scFactor = 0; % epsilon
    decoderN = 'svm';
    
    for s = 0:0.5:10
        test_acc_iter_wh = []; test_acc_iter_nh = []; test_acc_iter_wl = []; test_acc_iter_nl = [];
        Dir1 = []; Dir2 = []; Dir3 = []; Dir4 = [];

        for r = 1:100
            
            %%%%%%%%%%%%%% wide prior %%%%%%%%%%%%%%
            
            %for neuron-neuron correlation structure
            corr = [mean(Monkeys(1).NNcorr) mean(Monkeys(2).NNcorr)]; % high/low
                        
            if contrast == 'h'
                corr = corr(1);
                Fano = nanmean(Monkeys(1).Fanofactor);
                gGain = gamrnd(params1(1),params1(2),Npd,nMulti);
                sponAct = gamrnd(params2(1),params2(2),Npd,nMulti);
                dirSigma = gamrnd(params5(1),params5(2),Npd,nMulti);
                cparams1 = gamrnd(params7(1,1),params7(1,2),Npd,nMulti);
                cparams2 = gamrnd(params7(2,1),params7(2,2),Npd,nMulti);
                cparams3 = gamrnd(params7(3,1),params7(3,2),Npd,nMulti);
                cparams4 = NaN;
            elseif contrast == 'l'
                corr = corr(2);
                Fano = nanmean(Monkeys(2).Fanofactor);
                gGain = gamrnd(params3(1),params3(2),Npd,nMulti);
                sponAct = gamrnd(params4(1),params4(2),Npd,nMulti);
                dirSigma = gamrnd(params6(1),params6(2),Npd,nMulti);
                cparams1 = gamrnd(params8(1,1),params8(1,2),Npd,nMulti);
                cparams2 = gamrnd(params8(2,1),params8(2,2),Npd,nMulti);
                cparams3 = gamrnd(params8(3,1),params8(3,2),Npd,nMulti);
                cparams4 = NaN;
            end
            
            %--------------------- first direction ---------------------%
            stimTheta = 0; % stimulus motion direction
            [respMatrix1, param1, mResp1] = makeRespMatWithConstantNoise(Npd, nMulti, gGain, stimTheta, dirSigma, sponAct, numTrials, corr, Fano, cparams1, cparams2, cparams3, cparams4, decoderN);
            
            %--------------------- second direction ---------------------%
            stimTheta = s; % stimulus direction            
            [respMatrix2, param2, mResp2] = makeRespMatWithConstantNoise(Npd, nMulti, gGain, stimTheta, dirSigma, sponAct, numTrials, corr, Fano, cparams1, cparams2, cparams3, cparams4, decoderN);
            
            %--------------------- Preparing Dataset ---------------------%
            stimTheta = 0;
            x1 = respMatrix1;
            y1 = stimTheta;
            
            stimTheta = s;
            x2 = respMatrix2;
            if s == 0
                y2 = 1;
            else
                y2 = stimTheta;
            end
            
            %------------------- Binary Classification -------------------%
            X = [x1; x2];
            y = [repmat(y1,size(x1,1),1); repmat(y2,size(x1,1),1)];
            
            rand_num = randperm(size(X,1));
            X_train = X(rand_num(1:round(0.8*length(rand_num))),:);
            y_train = y(rand_num(1:round(0.8*length(rand_num))),:);
            
            X_test = X(rand_num(round(0.8*length(rand_num))+1:end),:);
            y_test = y(rand_num(round(0.8*length(rand_num))+1:end),:);
            
            %--------------------- Train & Test ---------------------%
            Md1 = fitclinear(X_train,y_train);
            test_accuracy = sum((predict(Md1,X_test) == y_test))/length(y_test)*100;
            if contrast == 'h'
                test_acc_iter_wh = [test_acc_iter_wh; test_accuracy];
            elseif contrast == 'l'
                test_acc_iter_wl = [test_acc_iter_wl; test_accuracy];
            end
            
            %%%%%%%%%%%%%% narrow prior %%%%%%%%%%%%%%
            
            %--------------------- first direction ---------------------%
            %Preparing Dataset
            stimTheta = 0;
            
            %Modulation matrix
            modulation = [];
            if contrast == 'h'
                regression = Monkeys(1).Regression(twin);
            elseif contrast == 'l'
                regression = Monkeys(2).Regression(twin);
            end
            
            for i = 1:length(param1.PD)
                PD = param1.PD(i);
                dirD = (stimTheta - PD)*pi/180.0;
                if dirD > pi
                    dirD = dirD - 2*pi;
                end
                if dirD < -pi
                    dirD = dirD + 2*pi;
                end
                m = regression;
                modulation = [modulation exp(m*abs(dirD*180/pi))]; %response reduction rate
            end
            
            x1 = respMatrix1.*modulation;
            y1 = stimTheta;
            
            %--------------------- second direction ---------------------%
            %Preparing Dataset
            stimTheta = s;
            
            %Modulation matrix
            modulation = [];
            for i = 1:length(param2.PD)
                PD = param2.PD(i);
                dirD = (stimTheta - PD)*pi/180.0;
                if dirD > pi
                    dirD = dirD - 2*pi;
                end
                if dirD < -pi
                    dirD = dirD + 2*pi;
                end
                m = regression;
                modulation = [modulation exp(m*abs(dirD*180/pi))]; %response reduction rate
            end
            
            x2 = respMatrix2.*modulation;
            if s == 0
                y2 = 1;
            else
                y2 = stimTheta;
            end

            %Binary Classification
            X = [x1; x2];
            y = [repmat(y1,size(x1,1),1); repmat(y2,size(x2,1),1)];
            
            rand_num = randperm(size(X,1));
            X_train = X(rand_num(1:round(0.8*length(rand_num))),:);
            y_train = y(rand_num(1:round(0.8*length(rand_num))),:);
            
            X_test = X(rand_num(round(0.8*length(rand_num))+1:end),:);
            y_test = y(rand_num(round(0.8*length(rand_num))+1:end),:);
            
            %Train & Test
            Md1 = fitclinear(X_train,y_train);
            test_accuracy = sum((predict(Md1,X_test) == y_test))/length(y_test)*100;
            if contrast == 'h'
                test_acc_iter_nh = [test_acc_iter_nh; test_accuracy];
            elseif contrast == 'l'
                test_acc_iter_nl = [test_acc_iter_nl; test_accuracy];
            end
        end
        
        if contrast == 'h'
            save(['d' num2str(s) '_w' contrast '_test_acc_withoutFS.mat'],'test_acc_iter_wh');
            save(['d' num2str(s) '_n' contrast '_test_acc_withoutFS.mat'],'test_acc_iter_nh');
        elseif contrast == 'l'
            save(['d' num2str(s) '_w' contrast '_test_acc_withoutFS.mat'],'test_acc_iter_wl');
            save(['d' num2str(s) '_n' contrast '_test_acc_withoutFS.mat'],'test_acc_iter_nl');
        end        
    end    
end

%% ttest
clear; clc;

contrast = 'h'; % h for high contrast, l for low contrast

for s = 0:0.5:10
    load(['d' num2str(s) '_w' contrast '_test_acc_withoutFS.mat']);
    load(['d' num2str(s) '_n' contrast '_test_acc_withoutFS.mat']);
       
    if contrast == 'h'
        test_acc_iter_w = test_acc_iter_wh;
        test_acc_iter_n = test_acc_iter_nh;
    elseif contrast == 'l'
        test_acc_iter_w = test_acc_iter_wl;
        test_acc_iter_n = test_acc_iter_nl;
    end

    ttestResult = [];
    [h,p,ci,stats] = ttest(test_acc_iter_w,test_acc_iter_n);
    ttestResult(1,1) = h;
    ttestResult(1,2) = p;
    ttestResult(1,3) = stats.tstat;
    ttestResult(1,4) = nanmean(test_acc_iter_w);
    ttestResult(1,5) = nanmean(test_acc_iter_n);
    
    save(['d' num2str(s) '_ttest_' contrast '.mat'],'ttestResult');
end

%% values for psychometric curve

contrast = 'h'; % 'h' for high contrast, 'l' for low contrast

mw = []; sw = []; mn = []; sn = []; ttest = [];
for degree = 0%:0.5:10

    load(['d' num2str(degree) '_w' contrast '_test_acc_withoutFS.mat']);
    load(['d' num2str(degree) '_n' contrast '_test_acc_withoutFS.mat']);
    load(['d' num2str(degree) '_ttest_' contrast '.mat']);
    
    if contrast == 'h'
        mw = [mw mean(test_acc_iter_wh)];
        sw = [sw std(test_acc_iter_wh)];
        mn = [mn mean(test_acc_iter_nh)];
        sn = [sn std(test_acc_iter_nh)];
        ttest = [ttest ttestResult(1,1)];
    elseif contrast == 'l'
        mw = [mw mean(test_acc_iter_wl)];
        sw = [sw std(test_acc_iter_wl)];
        mn = [mn mean(test_acc_iter_nl)];
        sn = [sn std(test_acc_iter_nl)];
        ttest = [ttest ttestResult(1,1)];
    end
end

save(['value_' contrast '_real_85.mat'],'mw','sw','mn','sn','ttest');

% cumulative gaussian fitting
x = 0:0.5:10;
y1 = mw/100;
y2 = mn/100;

if contrast == 'h'
test_acc_iter_w = test_acc_iter_wh;
test_acc_iter_n = test_acc_iter_nh;
else
test_acc_iter_w = test_acc_iter_wl;
test_acc_iter_n = test_acc_iter_nl;
end
s1 = sw./sqrt(length(test_acc_iter_w));
s2 = sn./sqrt(length(test_acc_iter_n));

m = mean(x);
sig = (x(end)-m)/3;

% cumulative gaussian function
f = @(p,x) (1/2 * (1+erf ((x - p(1)) ./ (p(2).*sqrt(2)) )))/2+0.5; %p1 = mean / p2 = std / /2+0.5 is for this data

% find p for minimum norm(=Euclidian length = sqrt(sum of squares))
P1 = fminsearch(@(p) norm(y1 - f(p,x)), [m sig]); 
P2 = fminsearch(@(p) norm(y2 - f(p,x)), [m sig]);

figure;
x_vector=min(x):(max(x)-min(x))/100:max(x);
plot(x,y1*100,'k.',x_vector,f(P1,x_vector)*100,'k-')
hold on;
plot(x,y2*100,'r.',x_vector,f(P2,x_vector)*100,'r-')
xlim([0 10]);

savefig(['cGaussian_' contrast '_85.fig'])
save(['paramsCgau_' contrast '_85.mat'],'P1','P2','y1','y2','s1','s2');

%% High and Low in a same figure

x = 0:0.5:10;
% cumulative gaussian function
f = @(p,x) (1/2 * (1+erf ((x - p(1)) ./ (p(2).*sqrt(2)) )))/2+0.5; %p1 = mean / p2 = std / /2+0.5 is for this data

fig = figure;
load('paramsCgau_h_85.mat');
x_vector=min(x):(max(x)-min(x))/100:max(x);
plot(x,y1*100,'b.',x_vector,f(P1,x_vector)*100,'b-')
hold on;
errorbar(x,y1*100,s1);
hold on;
plot(x,y2*100,'r.',x_vector,f(P2,x_vector)*100,'r-')
hold on;
errorbar(x,y2*100,s2);
hold on;
xlim([0 10]);
hold on;

load('paramsCgau_l_85.mat');
hold on;
plot(x,y1*100,'c.',x_vector,f(P1,x_vector)*100,'c:')
hold on;
errorbar(x,y1*100,s1);
hold on;
plot(x,y2*100,'m.',x_vector,f(P2,x_vector)*100,'m:')
hold on;
errorbar(x,y2*100,s2);
hold on;
xlim([0 10]);

savefig('cGaussian_85.fig')

%% params of psychometric curve
clear;

P1 = []; P2 = [];
contrast = 'h'; % 'h' for high contrast, 'l' for low contrast
for n = 1:100
    mw = []; mn = [];
    for degree = 0:0.5:10
        load(['d' num2str(degree) '_w' contrast '_test_acc_withoutFS.mat']);
        load(['d' num2str(degree) '_n' contrast '_test_acc_withoutFS.mat']);
        
        if contrast == 'h'
            mw = [mw test_acc_iter_wh(n)];
            mn = [mn test_acc_iter_nh(n)];
        elseif contrast == 'l'
            mw = [mw test_acc_iter_wl(n)];
            mn = [mn test_acc_iter_nl(n)];
        end
    end
    
    % cumulative gaussian fitting
    x = 0:0.5:10;
    y1 = mw/100; %percent to ratio
    y2 = mn/100;
    
    if contrast == 'h'
        test_acc_iter_w = test_acc_iter_wh;
        test_acc_iter_n = test_acc_iter_nh;
    elseif contrast == 'l'
        test_acc_iter_w = test_acc_iter_wl;
        test_acc_iter_n = test_acc_iter_nl;
    end
    
    m = mean(x);
    sig = (x(end)-m)/3;
    
    % cumulative gaussian function
    f = @(p,x) (1/2 * (1+erf ((x - p(1)) ./ (p(2).*sqrt(2)) )))/2+0.5; %p1 = mean / p2 = std / /2+0.5 is for this data
    
    %find p for minimum norm(=Euclidian length = sqrt(sum of squares))
    P1 = [P1; fminsearch(@(p) norm(y1 - f(p,x)), [m sig])]; %wide prior
    P2 = [P2; fminsearch(@(p) norm(y2 - f(p,x)), [m sig])]; %narrow prior
    
end

save(['params_' contrast '_85.mat'],'P1','P2');

%% ttest of params

% High contrast - wide vs narrow
clear;
load('params_h_85.mat');
ttestMean = [];
[h,p,ci,stats] = ttest(P1(:,1),P2(:,1));
ttestMean(1,1) = h;
ttestMean(1,2) = p;
ttestMean(1,3) = stats.tstat;
ttestMean(1,4) = nanmean(P1(:,1));
ttestMean(1,5) = nanmean(P2(:,1));

ttestStd = [];
[h,p,ci,stats] = ttest(P1(:,2),P2(:,2));
ttestStd(1,1) = h;
ttestStd(1,2) = p;
ttestStd(1,3) = stats.tstat;
ttestStd(1,4) = nanmean(P1(:,2));
ttestStd(1,5) = nanmean(P2(:,2));

ebP1m = std(P1(:,1))/sqrt(100);
ebP2m = std(P2(:,1))/sqrt(100);
ebP1s = std(P1(:,2))/sqrt(100);
ebP2s = std(P2(:,2))/sqrt(100);
save('params_ttest_h_85.mat','ttestMean','ttestStd','ebP1m','ebP2m','ebP1s','ebP2s');

fig = figure;
subplot(1,2,1)
bar(1:2,[ttestMean(1,4) ttestMean(1,5)]);
hold on;
errorbar(1:2,[ttestMean(1,4) ttestMean(1,5)],[ebP1m ebP2m]);
subplot(1,2,2)
bar(1:2,[ttestStd(1,4) ttestStd(1,5)]);
hold on;
errorbar(1:2,[ttestStd(1,4) ttestStd(1,5)],[ebP1s ebP2s]);
savefig('params_ttest_h_85.fig')

% Low contrast - wide vs narrow
clear;
load('params_l_85.mat');
ttestMean = [];
[h,p,ci,stats] = ttest(P1(:,1),P2(:,1));
ttestMean(1,1) = h;
ttestMean(1,2) = p;
ttestMean(1,3) = stats.tstat;
ttestMean(1,4) = nanmean(P1(:,1));
ttestMean(1,5) = nanmean(P2(:,1));

ttestStd = [];
[h,p,ci,stats] = ttest(P1(:,2),P2(:,2));
ttestStd(1,1) = h;
ttestStd(1,2) = p;
ttestStd(1,3) = stats.tstat;
ttestStd(1,4) = nanmean(P1(:,2));
ttestStd(1,5) = nanmean(P2(:,2));

ebP1m = std(P1(:,1))/sqrt(100);
ebP2m = std(P2(:,1))/sqrt(100);
ebP1s = std(P1(:,2))/sqrt(100);
ebP2s = std(P2(:,2))/sqrt(100);
save('params_ttest_L_85.mat','ttestMean','ttestStd','ebP1m','ebP2m','ebP1s','ebP2s');

figure;
subplot(1,2,1)
bar(1:2,[ttestMean(1,4) ttestMean(1,5)]);
hold on;
errorbar(1:2,[ttestMean(1,4) ttestMean(1,5)],[ebP1m ebP2m]);
subplot(1,2,2)
bar(1:2,[ttestStd(1,4) ttestStd(1,5)]);
hold on;
errorbar(1:2,[ttestStd(1,4) ttestStd(1,5)],[ebP1s ebP2s]);
savefig('params_ttest_L_85.fig')

% Comparison between high and low
clear;
load('params_h_85.mat');
hp1 = P1; hp2 = P2;
load('params_l_85.mat');
lp1 = P1; lp2 = P2;

%wide
ttestMean = [];
[h,p,ci,stats] = ttest(hp1(:,1),lp1(:,1));
ttestMean(1,1) = h;
ttestMean(1,2) = p;
ttestMean(1,3) = stats.tstat;
ttestMean(1,4) = nanmean(hp1(:,1));
ttestMean(1,5) = nanmean(lp1(:,1));

ttestStd = [];
[h,p,ci,stats] = ttest(hp1(:,2),lp1(:,2));
ttestStd(1,1) = h;
ttestStd(1,2) = p;
ttestStd(1,3) = stats.tstat;
ttestStd(1,4) = nanmean(hp1(:,2));
ttestStd(1,5) = nanmean(lp1(:,2));

ebP1m = std(P1(:,1))/sqrt(100);
ebP2m = std(P2(:,1))/sqrt(100);
ebP1s = std(P1(:,2))/sqrt(100);
ebP2s = std(P2(:,2))/sqrt(100);
save('params_ttest_HL_wide.mat','ttestMean','ttestStd','ebP1m','ebP2m','ebP1s','ebP2s');

figure;
subplot(1,2,1)
bar(1:2,[ttestMean(1,4) ttestMean(1,5)]);
hold on;
errorbar(1:2,[ttestMean(1,4) ttestMean(1,5)],[ebP1m ebP2m]);
subplot(1,2,2)
bar(1:2,[ttestStd(1,4) ttestStd(1,5)]);
hold on;
errorbar(1:2,[ttestStd(1,4) ttestStd(1,5)],[ebP1s ebP2s]);
savefig('params_ttest_HL_wide.fig')

%narrow
ttestMean = [];
[h,p,ci,stats] = ttest(hp2(:,1),lp2(:,1));
ttestMean(1,1) = h;
ttestMean(1,2) = p;
ttestMean(1,3) = stats.tstat;
ttestMean(1,4) = nanmean(hp2(:,1));
ttestMean(1,5) = nanmean(lp2(:,1));

ttestStd = [];
[h,p,ci,stats] = ttest(hp2(:,2),lp2(:,2));
ttestStd(1,1) = h;
ttestStd(1,2) = p;
ttestStd(1,3) = stats.tstat;
ttestStd(1,4) = nanmean(hp2(:,2));
ttestStd(1,5) = nanmean(lp2(:,2));

ebP1m = std(hp2(:,1))/sqrt(100);
ebP2m = std(lp2(:,1))/sqrt(100);
ebP1s = std(hp2(:,2))/sqrt(100);
ebP2s = std(lp2(:,2))/sqrt(100);
save('params_ttest_HL_narrow.mat','ttestMean','ttestStd','ebP1m','ebP2m','ebP1s','ebP2s');

figure;
subplot(1,2,1)
bar(1:2,[ttestMean(1,4) ttestMean(1,5)]);
hold on;
errorbar(1:2,[ttestMean(1,4) ttestMean(1,5)],[ebP1m ebP2m]);
subplot(1,2,2)
bar(1:2,[ttestStd(1,4) ttestStd(1,5)]);
hold on;
errorbar(1:2,[ttestStd(1,4) ttestStd(1,5)],[ebP1s ebP2s]);
savefig('params_ttest_HL_narrow.fig')