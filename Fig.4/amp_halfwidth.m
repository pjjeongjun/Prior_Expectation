%% Tuning amplitude and half-width
% Preferred direction was calculated using both wide and narrow prior block data.

load('AmpHW.mat');

% Tuning amplitude and half-width
whAmp = [data.A.whAmp; data.B.whAmp];
nhAmp = [data.A.nhAmp; data.B.nhAmp];
wlAmp = [data.A.wlAmp; data.B.wlAmp];
nlAmp = [data.A.nlAmp; data.B.nlAmp];
whHW = [data.A.whHW; data.B.whHW];
nhHW = [data.A.nhHW; data.B.nhHW];
wlHW = [data.A.wlHW; data.B.wlHW];
nlHW = [data.A.nlHW; data.B.nlHW];

ttestResultHAmp = [];
[h,p,ci,stats] = ttest(whAmp,nhAmp);
ttestResultHAmp(1,1) = h;
ttestResultHAmp(1,2) = p;
ttestResultHAmp(1,3) = stats.tstat;
ttestResultHAmp(1,4) = mean(whAmp);
ttestResultHAmp(1,5) = mean(nhAmp);

ttestResultLAmp = [];
[h,p,ci,stats] = ttest(wlAmp,nlAmp);
ttestResultLAmp(1,1) = h;
ttestResultLAmp(1,2) = p;
ttestResultLAmp(1,3) = stats.tstat;
ttestResultLAmp(1,4) = mean(wlAmp);
ttestResultLAmp(1,5) = mean(nlAmp);

ttestResultHHW = [];
[h,p,ci,stats] = ttest(whHW,nhHW);
ttestResultHHW(1,1) = h;
ttestResultHHW(1,2) = p;
ttestResultHHW(1,3) = stats.tstat;
ttestResultHHW(1,4) = mean(whHW);
ttestResultHHW(1,5) = mean(nhHW);

ttestResultLHW = [];
[h,p,ci,stats] = ttest(wlHW,nlHW);
ttestResultLHW(1,1) = h;
ttestResultLHW(1,2) = p;
ttestResultLHW(1,3) = stats.tstat;
ttestResultLHW(1,4) = mean(wlHW);
ttestResultLHW(1,5) = mean(nlHW);

save('AmpandHWP.mat','whAmp','nhAmp','wlAmp','nlAmp','ttestResultHAmp','ttestResultLAmp','whHW','nhHW','wlHW','nlHW','ttestResultHHW','ttestResultLHW');