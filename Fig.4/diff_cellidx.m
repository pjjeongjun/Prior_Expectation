%% larger and smaller delta-theta neurons
% indices for neurons based on the difference of preferred and prior directions
% Fig. 4, d-e

clear; clc;

load('pf_tar_dirs.mat');
pfDir = pfDir_A; % pfDir_A: monkey A, pfDir_B: monkey B
tarDir = tarDir_A; % tarDir_A: monkey A, tarDir_B: monkey B

percent = 0.3;

% difference between preferred and target directions
diff = [];
for i = 1:length(pfDir)
    if ~isnan(pfDir{i})
        if abs(pfDir{i}-tarDir{i}) > 180
            tmp = 360-abs(pfDir{i}-tarDir{i});
        else
            tmp = abs(pfDir{i}-tarDir{i});
        end
        diff = [diff; tmp];
    else
        diff= [diff; NaN];
    end
end

totalnum = length(find(~isnan(diff)));
num = round(totalnum*percent);
[ys,idxs] = sort(diff);
[yl,idxl] = sort(diff,'descend');
nanidx = max(find(isnan(yl)));

%Indices of cells
sidx = idxs(1:num); %small difference
lidx = idxl(nanidx+1:nanidx+1+num-1); %large difference
