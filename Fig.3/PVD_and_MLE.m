%% Population Vector Decoder / Maximum Likelihood Estimate
% Fig. 3, e-f

clear; clc; close all;

% load experimental data
load('Monkeys_cir.mat');
Monkeys = Data;

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
    twin = 18; %18 = 85ms
    nMulti = 10; % number of neurons for each direction
    Npd = 360; % number of prefer directions
    numTrials = 200; % number of trials (make it similar to acutual trial num)
    numSimul = 100;
    scFactor = 0*nMulti*Npd; % epsilon
    priorTheta = 0;
    decoderN = 'pvdmle';

    mle=struct('wide',zeros(numSimul,numTrials),'narrow',zeros(numSimul,numTrials));
    stimTheta = 0; % stimulus direction
    Dir1 = []; Dir2 = [];
    for r = 1:numSimul

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

        [respMatrix1, param, mResp1] = makeRespMatWithConstantNoise(Npd, nMulti, gGain, stimTheta, dirSigma, sponAct, numTrials, corr, Fano, cparams1, cparams2, cparams3, cparams4, decoderN);
        % PVD
        tmp1 = simpleDirVA(respMatrix1, param, scFactor); % decode direction
        Dir1 = [Dir1 tmp1];

        % MLE
        range=-179:1:180;
        pds=deg2rad(param.PD);
        aLL=[];
        for deg=range
            d1=deg2rad(0); d2=deg2rad(45); d3=deg2rad(deg); dLL=[];
            for tr=1:size(respMatrix1,1)
                LL=[cos(d1) sin(d1); cos(d2) sin(d2)]*[sum(respMatrix1(tr,:).*cos(pds)); sum(respMatrix1(tr,:).*sin(pds))];
                if deg==0
                    nLL=LL(1);
                elseif deg==45
                    nLL=LL(2);
                else
                    nLL=[cos(d3) sin(d3)]*inv([cos(d1) sin(d1); cos(d2) sin(d2)])*[LL(1); LL(2)];
                end
                dLL=[dLL nLL];
            end
            aLL=[aLL; dLL]; %deg*trials
        end

        mle_w=[];
        for tr=1:size(aLL,2)
            [m,d]=max(aLL(:,tr));
            mle_w=[mle_w range(d)];
        end

        %%%%%%%%%%%%%% narrow prior %%%%%%%%%%%%%%

        %Modulation matrix
        modulation = [];
        if contrast == 'h'
            regression = Monkeys(1).Regression(twin); %%%%% high
        elseif contrast == 'l'
            regression = Monkeys(2).Regression(twin);
        end

        for i = 1:length(param.PD)
            PD = param.PD(i);
            dirD = (priorTheta - PD)*pi/180.0;
            if dirD > pi
                dirD = dirD - 2*pi;
            end
            if dirD < -pi
                dirD = dirD + 2*pi;
            end
            modulation = [modulation exp(regression*abs(dirD*180/pi))]; %response reduction rate
        end

        respMatrix2 = respMatrix1.*modulation;
        % PVD
        tmp2 = simpleDirVA(respMatrix2, param, scFactor); % decode direction
        Dir2 = [Dir2 tmp2];

        % MLE
        aLL=[];
        for deg=range
            d1=deg2rad(0); d2=deg2rad(45); d3=deg2rad(deg); dLL=[];
            for tr=1:size(respMatrix2,1)
                LL=[cos(d1) sin(d1); cos(d2) sin(d2)]*[sum(respMatrix2(tr,:).*cos(pds)); sum(respMatrix2(tr,:).*sin(pds))];
                if deg==0
                    nLL=LL(1);
                elseif deg==45
                    nLL=LL(2);
                else
                    nLL=[cos(d3) sin(d3)]*inv([cos(d1) sin(d1); cos(d2) sin(d2)])*[LL(1); LL(2)];
                end
                dLL=[dLL nLL];
            end
            aLL=[aLL; dLL]; %deg*trials
        end

        mle_n=[];
        for tr=1:size(aLL,2)
            [m,d]=max(aLL(:,tr));
            mle_n=[mle_n range(d)];
        end
        mle.wide(r,:)=mle_w;
        mle.narrow(r,:)=mle_n;
    end

    %0 degree direction
    decoder.mle=mle;
    decoder.pvd.wide=Dir1;
    decoder.pvd.narrow=Dir2;
    save(['' contrast '_decoder.mat'],'decoder');
end

%% SD ratio from simulation

stimTheta = 0;
contrast = 'h'; % h for high, l fo low
decoderN = 'PVD'; % PVD or MLE

load(['' contrast '_decoder.mat']);

m1=[];m2=[];s1=[];s2=[];

for r = 1:size(decoder.mle.wide,1)
    %mean
    if strcmp(decoderN,'PVD')
        m1 = [m1; mean(decoder.pvd.wide(r,:))];
        m2 = [m2; mean(decoder.pvd.narrow(r,:))];
    elseif strcmp(decoderN,'MLE')
        m1 = [m1; mean(decoder.mle.wide(r,:))];
        m2 = [m2; mean(decoder.mle.narrow(r,:))];
    end
    %std
    if strcmp(decoderN,'PVD')
        s1 = [s1; std(decoder.pvd.wide(r,:))];
        s2 = [s2; std(decoder.pvd.narrow(r,:))];
    elseif strcmp(decoderN,'MLE')
        s1 = [s1; std(decoder.mle.wide(r,:))];
        s2 = [s2; std(decoder.mle.narrow(r,:))];
    end
end

%std ratio
sr = s2./s1;

%ttest
ttestResult = [];
[h,p,ci,stats] = ttest(sr);
ttestResult(1,1) = h;
ttestResult(1,2) = p;
ttestResult(1,3) = stats.tstat;
ttestResult(1,4) = mean(sr);

save(['d' num2str(stimTheta) '_ttest_stdratio_' contrast '_' decoderN '.mat'],'sr','ttestResult');

%% Comparison of SD ratio between experimental and simulation
clear; clc;

decoderN = 'PVD'; %PVD or MLE

load('stdratio.mat');
if strcmp(decoderN,'PVD')
    load('d0_ttest_stdratio_l_pvd.mat');
    lowsr = sr;
    load('d0_ttest_stdratio_h_pvd.mat');
    highsr = sr;
elseif strcmp(decoderN,'MLE')
    load('d0_ttest_stdratio_l_mle.mat');
    lowsr = sr;
    load('d0_ttest_stdratio_h_mle.mat');
    highsr = sr;
end

fig = figure;
scatter(high,low) %actual
hold on;
scatter(highsr,lowsr) %simulated
hold on;
errorbar(mean(high),mean(low),std(high),'horizontal','blue')
hold on;
errorbar(mean(high),mean(low),std(low),'vertical','blue')
hold on;
errorbar(mean(highsr),mean(lowsr),std(highsr),'horizontal','red')
hold on;
errorbar(mean(highsr),mean(lowsr),std(lowsr),'vertical','red')
savefig(['' decoderN '.fig'])

%ttestH - std log ratio
ttestResultH = [];
[h,p,ci,stats] = ttest(highsr);
ttestResultH(1,1) = h;
ttestResultH(1,2) = p;
ttestResultH(1,3) = stats.tstat;
ttestResultH(1,4) = mean(highsr);

%ttestL - std log ratio
ttestResultL = [];
[h,p,ci,stats] = ttest(lowsr);
ttestResultL(1,1) = h;
ttestResultL(1,2) = p;
ttestResultL(1,3) = stats.tstat;
ttestResultL(1,4) = mean(lowsr);

save(['stdratio_simulated_' decoderN '.mat'],'highsr','lowsr','ttestResultH','ttestResultL');

%ttestH - std log ratio
ttestResultH = [];
[h,p,ci,stats] = ttest2(high,highsr);
ttestResultH(1,1) = h;
ttestResultH(1,2) = p;
ttestResultH(1,3) = stats.tstat;
ttestResultH(1,4) = mean(high);
ttestResultH(1,5) = mean(highsr);

%ttestL - std log ratio
ttestResultL = [];
[h,p,ci,stats] = ttest2(low,lowsr);
ttestResultL(1,1) = h;
ttestResultL(1,2) = p;
ttestResultL(1,3) = stats.tstat;
ttestResultL(1,4) = mean(low);
ttestResultL(1,5) = mean(lowsr);

%ttest - std log ratio
ttestResult = [];
[h,p,ci,stats] = ttest(highsr,lowsr);
ttestResult(1,1) = h;
ttestResult(1,2) = p;
ttestResult(1,3) = stats.tstat;
ttestResult(1,4) = mean(highsr);
ttestResult(1,5) = mean(lowsr);

save(['stdratio_comparison_' decoderN '.mat'],'ttestResultH','ttestResultL','ttestResult');
