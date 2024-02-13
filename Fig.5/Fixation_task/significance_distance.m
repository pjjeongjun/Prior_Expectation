%% significance of distance
clear; clc;

%%%% distance in subspace %%%

%load data
load('subspace_position_A.mat');
real_a = sdata;
load('subspace_position_B.mat');
real_b = sdata;
load('permutation_A.mat');
perm_a = sdata;
load('permutation_B.mat');
perm_b = sdata;

% null
distance = struct('A',[],'B',[]);
distance.A = struct('high',[],'low',[]);
distance.B = struct('high',[],'low',[]);
for r = 1:1000 %samples
    ah_dis = []; al_dis = []; bh_dis = []; bl_dis = [];
    for t = 1:27 %time=101:20:621;(window size:20ms)
        %High
        %A
        wx = perm_a.high.wide.x(:,t,r); wy = perm_a.high.wide.y(:,t,r);
        nx = perm_a.high.narrow.x(:,t,r); ny = perm_a.high.narrow.y(:,t,r);
        ah_dis = [ah_dis nanmean(wx-nx)];
        %B
        wx = perm_b.high.wide.x(:,t,r); wy = perm_b.high.wide.y(:,t,r);
        nx = perm_b.high.narrow.x(:,t,r); ny = perm_b.high.narrow.y(:,t,r);
        bh_dis = [bh_dis nanmean(wx-nx)];

        %Low
        %A
        wx = perm_a.low.wide.x(:,t,r); wy = perm_a.low.wide.y(:,t,r);
        nx = perm_a.low.narrow.x(:,t,r); ny = perm_a.low.narrow.y(:,t,r);
        al_dis = [al_dis nanmean(wx-nx)];

        %B
        wx = perm_b.low.wide.x(:,t,r); wy = perm_b.low.wide.y(:,t,r);
        nx = perm_b.low.narrow.x(:,t,r); ny = perm_b.low.narrow.y(:,t,r);
        bl_dis = [bl_dis nanmean(wx-nx)];
    end
    distance.A.high = [distance.A.high; ah_dis];
    distance.A.low = [distance.A.low; al_dis];
    distance.B.high = [distance.B.high; bh_dis];
    distance.B.low = [distance.B.low; bl_dis];
end
nulldist = distance;

% real
load('distance_data.mat');
realdist = distance;

% permutation test for specific duration
% period = 6:10; % -100~0
period = 6:20; % -100~300

% %one-tailed
% real = mean(realdist.A.high(period));
% null = mean(nulldist.A.high(:,period),2);
% ah_p = (1+length(find(null>=real)))/(1+1000)
% 
% real = mean(realdist.A.low(period));
% null = mean(nulldist.A.low(:,period),2);
% al_p = (1+length(find(null>=real)))/(1+1000)
% 
% real = mean(realdist.B.high(period));
% null = mean(nulldist.B.high(:,period),2);
% bh_p = (1+length(find(null>=real)))/(1+1000)
% 
% real = mean(realdist.B.low(period));
% null = mean(nulldist.B.low(:,period),2);
% bl_p = (1+length(find(null>=real)))/(1+1000)

%two-tailed
ah = mean(realdist.A.high(period))
null = mean(nulldist.A.high(:,period),2);
real = ah-mean(null);
null = null-mean(null);
ah_p = (1+length(find(abs(null)>=abs(real))))/(1+length(null))

al = mean(realdist.A.low(period))
null = mean(nulldist.A.low(:,period),2);
real = al-mean(null);
null = null-mean(null);
al_p = (1+length(find(abs(null)>=abs(real))))/(1+length(null))

bh = mean(realdist.B.high(period))
null = mean(nulldist.B.high(:,period),2);
real = bh-mean(null);
null = null-mean(null);
bh_p = (1+length(find(abs(null)>=abs(real))))/(1+length(null))

bl = mean(realdist.B.low(period))
null = mean(nulldist.B.low(:,period),2);
real = bl-mean(null);
null = null-mean(null);
bl_p = (1+length(find(abs(null)>=abs(real))))/(1+length(null))
