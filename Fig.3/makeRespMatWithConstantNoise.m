function [respMatrix, params, mResp] = makeRespMatWithConstantNoise(Npd, numRow, gain, direction, sigD, sponR, nTrials, corr, Fano, cparams1, cparams2, cparams3, cparams4, decoderN)
 
[mResp, params, paramsRaw] = makeMeanResp(Npd, numRow, gain, direction, sigD, sponR, cparams1, cparams2, cparams3, cparams4, decoderN);

nParams.PD = (params.PD)*pi./180.0;
nParams.Multi = params.Multi;
nParams.nRows = numRow;
nParams.tDir = direction;

n = length(mResp.vec);

% correlation matrix
corrMat = ones(Npd*numRow,Npd*numRow).*corr; % select from data-driven distribution
for i = 1:n
	corrMat(i,i) = 1;
end

numTrials = nTrials;
yMatrix = zeros(numTrials, n);

ranMatrix = normrnd(0, 1, numTrials, n);

Q = chol(corrMat, 'lower'); % Do the cholesky decomposition for getting square root of correlation matrix
for i = 1:numTrials
	yMatrix(i,:) = ranMatrix(i,:)*transpose(Q);
end

correctionFactor = mean(var(ranMatrix))/mean(var(yMatrix));

meanYMatrix = mean(yMatrix);

respMatrix = zeros(numTrials, n);

for i = 1:numTrials
	respMatrix(i,:) = mResp.vec + (yMatrix(i, :)-meanYMatrix).*sqrt(mResp.vec.*correctionFactor.*Fano);
end