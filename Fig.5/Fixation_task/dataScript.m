%% Reconstruct data for all conditions
clear; clc;

t = 101:20:621;
evCrit = 0.5;

for m = 1:2
    if m == 1
        load(['Data_A_trial_fixation_ev' num2str(evCrit) '.mat']);
    else
        load(['Data_B_trial_fixation_ev' num2str(evCrit) '.mat']);
    end
    % firing rate

    dataT = struct('unit',[],'time',[]);

    %unit
    dataT.unit = struct('response',[],'task_variable',[],'dimension',[]);
    for i = 1:size(c1,1) %neuron
        tmp1 = []; tmp2 = []; tmp3 = []; tmp4 = [];
        for d = 1:12 %direction
            dir = 30*(d-1);
            % responses in each condition
            tmp1 = [tmp1; c1{i,d}; c2{i,d}; c3{i,d}; c4{i,d}];

            % target dir in radians
            tmp2 = [tmp2; repmat(deg2rad(dir),[size(c1{i,d},1),1]); repmat(deg2rad(dir),[size(c2{i,d},1),1]); ...
                repmat(deg2rad(dir),[size(c3{i,d},1),1]); repmat(deg2rad(dir),[size(c4{i,d},1),1])];

            % 1 = wide prior / 0 = narrow prior
            tmp3 = [tmp3; repmat(1,[size(c1{i,d},1),1]); repmat(0,[size(c2{i,d},1),1]); ...
                repmat(1,[size(c3{i,d},1),1]); repmat(0,[size(c4{i,d},1),1])];

            % 1 = high contrast / 0 = low contrast
            tmp4 = [tmp4; repmat(1,[size(c1{i,d},1),1]); repmat(1,[size(c2{i,d},1),1]); ...
                repmat(0,[size(c3{i,d},1),1]); repmat(0,[size(c4{i,d},1),1]);];

            % unit number
            if i < 10
                dataT.unit(i).dimension = ['unit_00' num2str(i)''];
            elseif i < 100
                dataT.unit(i).dimension = ['unit_0' num2str(i)''];
            else
                dataT.unit(i).dimension = ['unit_' num2str(i)''];
            end
        end
        dataT.unit(i).response = tmp1;
        dataT.unit(i).task_variable.stim_dir = tmp2;
        dataT.unit(i).task_variable.prior = tmp3;
        dataT.unit(i).task_variable.contrast = tmp4;     
    end

    didx = [];
    for n = 1:size(dataT.unit,2)
        if isempty(dataT.unit(n).response)
            didx = [didx; n];
        end
    end

    dataT.unit(didx) = []; 
    dataT.unit = dataT.unit';

    %time
    dataT.time = t;
    if m == 1
        save('data_A.mat','dataT');
    else
        save('data_B.mat','dataT');
    end
end
