%% distance in subspace
clear; clc;

%load data

load('subspace_position_A.mat');
a = sdata;
load('subspace_position_B.mat');
b = sdata;

%
ah_dis = []; al_dis = []; bh_dis = []; bl_dis = [];
ah_wx = []; ah_nx = []; al_wx = []; al_nx = [];
bh_wx = []; bh_nx = []; bl_wx = []; bl_nx = [];

for n = 1:27 %time=101:20:621(window size:20ms)
    %High
    %A
    wx = a.high.wide.x(:,n); wy = a.high.wide.y(:,n);
    nx = a.high.narrow.x(:,n); ny = a.high.narrow.y(:,n);
    
    ah_dis = [ah_dis mean(wx-nx)];
    ah_wx = [ah_wx wx];
    ah_nx = [ah_nx nx];

    %B
    wx = b.high.wide.x(:,n); wy = b.high.wide.y(:,n);
    nx = b.high.narrow.x(:,n); ny = b.high.narrow.y(:,n);
    bh_dis = [bh_dis mean(wx-nx)];
    bh_wx = [bh_wx wx];
    bh_nx = [bh_nx nx];

    %Low
    %A
    wx = a.low.wide.x(:,n); wy = a.low.wide.y(:,n);
    nx = a.low.narrow.x(:,n); ny = a.low.narrow.y(:,n);
    al_dis = [al_dis mean(wx-nx)];
    al_wx = [al_wx wx];
    al_nx = [al_nx nx];

    %B
    wx = b.low.wide.x(:,n); wy = b.low.wide.y(:,n);
    nx = b.low.narrow.x(:,n); ny = b.low.narrow.y(:,n);
    bl_dis = [bl_dis mean(wx-nx)];
    bl_wx = [bl_wx wx];
    bl_nx = [bl_nx nx];
end

distance.A.high = ah_dis;
distance.A.low = al_dis;
distance.B.high = bh_dis;
distance.B.low = bl_dis;

save('distance_data.mat','distance','ah_wx','ah_nx','al_wx','al_nx','bh_wx','bh_nx','bl_wx','bl_nx');

t = 101:20:621;
t = t-300;
period = 1:26;

% figure
figure;
plot(t(period),ah_dis(period),'k');
hold on;
plot(t(period),al_dis(period),'r');
hold on;
plot(t(period),zeros(1,length(period)),'b.');
legend('high','low')
xlabel('time(ms)');
title('Distance between the two priors, A')
savefig('distance_A.fig');

figure;
plot(t(period),bh_dis(period),'k');
hold on;
plot(t(period),bl_dis(period),'r');
hold on;
plot(t(period),zeros(1,length(period)),'b.');
legend('high','low')
xlabel('time(ms)'); ylabel('distance');
title('Distance between the two priors, B')
savefig('distance_B.fig');

%% ttest
clear; clc;

load('distance_data.mat');

d1 = ah_wx-ah_nx;
d2 = al_wx-al_nx;
d3 = bh_wx-bh_nx;
d4 = bl_wx-bl_nx;

dd1 = []; dd2 = []; dd3 = []; dd4 = [];
for tt = 1:26
    dd1 = [dd1; d1(:,tt)];
    dd2 = [dd2; d2(:,tt)];
    dd3 = [dd3; d3(:,tt)];
    dd4 = [dd4; d4(:,tt)];
end

[h1,p1,ci1,stats1] = ttest(dd1);
[h2,p2,ci2,stats2] = ttest(dd2);
[h3,p3,ci3,stats3] = ttest(dd3);
[h4,p4,ci4,stats4] = ttest(dd4);

ttestresult = [];
ttestresult(:,1) = [h1; h2; h3; h4];
ttestresult(:,2) = [p1; p2; p3; p4];
ttestresult(:,3) = [stats1.tstat; stats2.tstat; stats3.tstat; stats4.tstat];
ttestresult(:,4) = [mean(dd1);mean(dd2);mean(dd3);mean(dd4)];

save('fixation.mat','dd1','dd2','dd3','dd4');