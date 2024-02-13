function [modelMeanResp, modelParamsVec, modelParams] = makeMeanResp(Npd, numRow, gain, direction, sigD, sponR, cparams1, cparams2, cparams3, cparams4, decoderN)

[modelParams] = makeParamMat(Npd, numRow, direction, decoderN);

modelMeanResp = [];
modelParamsVec.PD = zeros(1, length(modelParams.PD)*length(modelParams.Multi)); %prefer direction(PD) of the neuron
modelParamsVec.Multi = zeros(1, length(modelParams.PD)*length(modelParams.Multi)); %index of the neuron
modelMeanResp.vec = zeros(1, length(modelParams.PD)*length(modelParams.Multi)); %response of the neuron

k = 1;
for i = 1:length(modelParams.PD)
    for j = 1:length(modelParams.Multi)
        modelParamsVec.PD(k) = modelParams.PD(i);
        modelParamsVec.Multi(k) = modelParams.Multi(j);
        %Gaussian
        if isnan(cparams1(1)) && isnan(cparams2(1)) && isnan(cparams3(1)) && isnan(cparams4(1)) 
            modelMeanResp.vec(k) = mRespIndividual(gain(i,j), direction, modelParams.PD(i), sigD(i,j), sponR(i,j), cparams1, cparams2, cparams3);
            k = k+1;
        elseif ~isnan(cparams1(1)) && ~isnan(cparams2(1)) && ~isnan(cparams3(1)) && ~isnan(cparams4(1))
            modelMeanResp.vec(k) = mRespIndividual(gain(i,j), direction, modelParams.PD(i), sigD(i,j), sponR(i,j), cparams1(i,j), cparams2(i,j), cparams3(i,j), cparams4(i,j));
            k = k+1;
        %Circular 
        elseif ~isnan(cparams1(1)) && ~isnan(cparams2(1)) && ~isnan(cparams3(1)) && isnan(cparams4(1)) 
            modelMeanResp.vec(k) = mRespIndividual(gain(i,j), direction, modelParams.PD(i), sigD(i,j), sponR(i,j), cparams1(i,j), cparams2(i,j), cparams3(i,j), cparams4);
            k = k+1;
        end
    end
end
