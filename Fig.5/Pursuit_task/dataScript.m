%% Reshape data for all directions (prior direction, two outer directions)
clear; clc;

for m = 1:2
    if m == 1
        load('Data_A_trial.mat');
    else
        load('Data_B_trial.mat');
    end

    dataT = struct('unit',[],'time',[]);

    %unit
    dataT.unit = struct('response',[],'task_variable',[],'dimension',[]);
    for i = 1:size(c1,2)
        % responses in each condition
        dataT.unit(i).response = [c1{1,i}; c2{1,i}; c3{1,i}; c4{1,i}; ...
            c5{1,i}; c6{1,i}; c7{1,i}; c8{1,i}; c9{1,i}; c10{1,i}; c11{1,i}; c12{1,i}];

        % target dir in radians
        dataT.unit(i).task_variable.stim_dir = [repmat(deg2rad(0),[size(c1{1,i},1),1]); repmat(deg2rad(0),[size(c2{1,i},1),1]); ...
            repmat(deg2rad(0),[size(c3{1,i},1),1]); repmat(deg2rad(0),[size(c4{1,i},1),1]); repmat(deg2rad(-120),[size(c5{1,i},1),1]);...
            repmat(deg2rad(-120),[size(c6{1,i},1),1]); repmat(deg2rad(120),[size(c7{1,i},1),1]); repmat(deg2rad(120),[size(c8{1,i},1),1]); ...
            repmat(deg2rad(-15),[size(c9{1,i},1),1]); repmat(deg2rad(-15),[size(c10{1,i},1),1]); repmat(deg2rad(15),[size(c11{1,i},1),1]); repmat(deg2rad(15),[size(c12{1,i},1),1])];

        % 0 = wide prior / 1 = narrow prior
        dataT.unit(i).task_variable.prior = [repmat(0,[size(c1{1,i},1),1]); repmat(1,[size(c2{1,i},1),1]); ...
            repmat(0,[size(c3{1,i},1),1]); repmat(1,[size(c4{1,i},1),1]); repmat(0,[size(c5{1,i},1),1]);...
            repmat(0,[size(c6{1,i},1),1]); repmat(0,[size(c7{1,i},1),1]); repmat(0,[size(c8{1,i},1),1]); ...
            repmat(1,[size(c9{1,i},1),1]); repmat(1,[size(c10{1,i},1),1]); repmat(1,[size(c11{1,i},1),1]); repmat(1,[size(c12{1,i},1),1])];

        % 1 = high contrast / 0 = low contrast
        dataT.unit(i).task_variable.contrast = [repmat(1,[size(c1{1,i},1),1]); repmat(1,[size(c2{1,i},1),1]); ...
            repmat(0,[size(c3{1,i},1),1]); repmat(0,[size(c4{1,i},1),1]); repmat(1,[size(c5{1,i},1),1]);...
            repmat(0,[size(c6{1,i},1),1]); repmat(1,[size(c7{1,i},1),1]); repmat(0,[size(c8{1,i},1),1]); ...
            repmat(1,[size(c9{1,i},1),1]); repmat(0,[size(c10{1,i},1),1]); repmat(1,[size(c11{1,i},1),1]); repmat(0,[size(c12{1,i},1),1])];

        % unit number
        if i < 10
            dataT.unit(i).dimension = ['unit_00' num2str(i)''];
        elseif i < 100
            dataT.unit(i).dimension = ['unit_0' num2str(i)''];
        else
            dataT.unit(i).dimension = ['unit_' num2str(i)''];
        end
    end
    dataT.unit = dataT.unit';

    %time
    dataT.time = 1:20:981;
    if m == 1
        save('data_A.mat','dataT');
    else
        save('data_B.mat','dataT');
    end
end