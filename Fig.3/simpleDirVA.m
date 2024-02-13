function [Dir] = simpleDirVA(respMatrix, params, epsilon)

numTrials = size(respMatrix, 1);
numCells = size(respMatrix, 2);

pDir = params.PD;
Dir = NaN(numTrials, 1);
for i = 1:numTrials
	Dir(i) = sum(pDir.*respMatrix(i,:))./(epsilon+sum(respMatrix(i,:)));
end

%Dir = Dir+normrnd(0,1, numTrials,1);