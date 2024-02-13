function [modelParams] = makeParamMat(Npd, numDuple, direction, decoderN)

if strcmp(decoderN,'pvd')
    if direction==0
        tmpVal = linspace(-180, 180, Npd+1);
        modelParams.PD = tmpVal(2:end);
    elseif direction==-15
        modelParams.PD = linspace(-195, 165, Npd);
    elseif direction==15
        modelParams.PD = linspace(-165, 195, Npd);
    end
else
    tmpVal = linspace(-180, 180, Npd+1);
    modelParams.PD = tmpVal(2:end);
end

modelParams.Multi = 1:numDuple;
