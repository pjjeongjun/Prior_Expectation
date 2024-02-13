function dirTuningInfo = dirTuningEst(directions, resp, Wgt);

	rMax = max(resp);

	maxIndex = find(resp == rMax);
	if length(maxIndex) > 1
		maxIndex = maxIndex(1);
	end
	startDir = directions(maxIndex)*pi/180.0;
	startDirPre = directions(maxIndex)-60;
	if startDirPre < 0
		startDirPre = startDirPre+360.0;
	end
	startDirPre = startDirPre*pi/180.0;
	startDirPost = directions(maxIndex)+60;
	if startDirPost > 360.0
		startDirPost = startDirPost-360.0;
	end
	startDirPost = startDirPost*pi/180.0;
	
	x0 = [rMax 0.5 startDir;rMax 1 startDir;rMax 4 startDir;rMax 0.5 startDirPre;rMax 0.5 startDirPost];

	xData = directions*pi./180.0;
	yData = resp;
	if exist('Wgt')
		weight = Wgt;
		wZero = [];
		wZero = find(weight == 0);
		if ~isempty(wZero)
			weight(wZero) = 1;
		end
	else
		weight = [];
	end
	
	var0 = var(yData);
	x = zeros(size(x0));
	output = zeros(size(x0, 1), 1);
	
	for i = 1:size(x0, 1)
		if ~isempty(weight)
			[x(i,:),fval,eFlag,output(i)] = fitSumDirTuning(x0(i,:), xData, yData, weight);
		else
			[x(i,:),fval,eFlag,output(i)] = fitSumDirTuning(x0(i,:), xData, yData);
		end
	end
	maxOut = nanmax(output);
	maxOutArray = find(output == maxOut);
	ctDirs = [0:0.1:360];
	ctDirsRad = ctDirs*pi./180.0;
	if isempty(maxOutArray)
		param0 = [nan nan nan];
		ev0 = nan;
		estDir = param0;
		evDir = ev0;		
		tuningEstimates = NaN(1, length(ctDirs));
		offset = nan;
		peak = nan;
		halfWidth = nan;
	else
		param0 = x(maxOutArray(1), :);
		ev0 = maxOut;
		estDir = param0;
		evDir = ev0;		
		tuningEstimates = circularNormalFunc(estDir, ctDirsRad);
		offset = min(tuningEstimates);		
		peak = max(tuningEstimates) - offset;
		halfWidth = 2*acos((log(0.5 + 0.5*exp(-2*estDir(2))))/estDir(2)+1);
		halfWidth = halfWidth*180.0/pi;
	end

	dirTuningInfo = [];
	dirTuningInfo.params.circularGaussian = estDir;
	dirTuningInfo.params.ev = evDir;
	dirTuningInfo.estimates.offset = offset;
	dirTuningInfo.estimates.amplitude = peak;
	dirTuningInfo.estimates.halfWidth = halfWidth;
	dirTuningInfo.estimates.pfDir = estDir(3)*180.0/pi;
	dirTuningInfo.estimates.directions = ctDirs;
	dirTuningInfo.estimates.tuningCurve = tuningEstimates;
	
	rMax = max(resp) - min(resp);
	maxIndex = find(resp == max(resp));
	if length(maxIndex) > 1
		maxIndex = maxIndex(1);
	end
	startDir = directions(maxIndex)*pi/180.0;
	startDirPre = directions(maxIndex)-60;
	if startDirPre < 0
		startDirPre = startDirPre+360.0;
	end
	startDirPre = startDirPre*pi/180.0;
	startDirPost = directions(maxIndex)+60;
	if startDirPost > 360.0
		startDirPost = startDirPost-360.0;
	end
	startDirPost = startDirPost*pi/180.0;
	
	x0 = [rMax startDir pi/4 min(resp);rMax startDir pi/2 min(resp);rMax startDir pi min(resp);rMax startDirPre pi/4 min(resp);rMax startDirPost pi/4 min(resp)];

	x = zeros(size(x0));
	output = zeros(size(x0, 1), 1);
	
	for i = 1:size(x0, 1)
		if ~isempty(weight)
			[x(i,:),fval,eFlag,output(i)] = fitSumDirGaussian(x0(i,:), xData, yData, weight);
		else
			[x(i,:),fval,eFlag,output(i)] = fitSumDirGaussian(x0(i,:), xData, yData);
		end
	end
	maxOut = nanmax(output);
	maxOutArray = find(output == maxOut);
	ctDirs = [0:0.1:360];
	ctDirsRad = ctDirs*pi./180.0;
	if isempty(maxOutArray)
		param0 = [nan nan nan nan];
		ev0 = nan;
		estDir = param0;
		evDir = ev0;		
		tuningEstimates = NaN(1, length(ctDirs));
		offset = nan;
		peak = nan;
		halfWidth = nan;
	else
		param0 = x(maxOutArray(1), :);
		ev0 = maxOut;
		estDir = param0;
		evDir = ev0;		
		tuningEstimates = gaussianTunning(estDir, ctDirsRad);
		offset = estDir(4);
		peak = estDir(1);
		halfWidth = 2*sqrt(2*log(2))*estDir(3);
		halfWidth = halfWidth*180.0/pi;
	end
	
	dirTuningInfo.params.gaussian.est = estDir;
	dirTuningInfo.params.gaussian.ev = evDir;
	dirTuningInfo.estimates.gaussian.offset = offset;
	dirTuningInfo.estimates.gaussian.amplitude = peak;
	dirTuningInfo.estimates.gaussian.halfWidth = halfWidth;
	dirTuningInfo.estimates.gaussian.pfDir = estDir(2)*180.0/pi;
	dirTuningInfo.estimates.gaussian.directions = ctDirs;
	dirTuningInfo.estimates.gaussian.tuningCurve = tuningEstimates;
	
function [x,fval,eFlag,output] = fitSumDirTuning(initGuess, directions, y, weight)

	x0 = initGuess;
	
	vlb = [1 0.001 0];                	% lower bounds
	vub = [800 30 2*pi];         		% upper bounds

	opts = optimset('fmincon');
	opts = optimset(opts, 'MaxFunEvals', 50000, 'MaxIter', 50000, 'TolFun', 1e-08, 'TolX', 1e-08, 'FunValCheck','on');
	opts = optimset(opts, 'LargeScale', 'off', 'Display', 'off', 'Algorithm', 'active-set');

	try
		if ~exist('weight')
			[x,fval,eFlag,output] = fmincon(@(param) sumLsqDirTuning(param, directions, y), x0, [], [], [], [], vlb, vub, [], opts);
		elseif exist('weight')
			[x,fval,eFlag,output] = fmincon(@(param) sumLsqDirTuning(param, directions, y, weight), x0, [], [], [], [], vlb, vub, [], opts);
		end
			
		y_hat = circularNormalFunc(x, directions);

		ss0 = sumLsqDirTuning(x, directions, y);
		ss0 = ss0/(length(y)-1);
		var0 = var(y);
		evVal = 1-ss0/var0;
		output = evVal;
	catch
		try
			opts1 = optimset('fminsearch');
			opts2 = optimset('fmincon');
			opts1 = optimset(opts1, 'MaxFunEvals', 50000, 'MaxIter', 50000, 'TolFun', 1e-012, 'TolX', 1e-012, 'FunValCheck','on');
			opts1 = optimset(opts1, 'LargeScale', 'off', 'Display', 'off', 'Algorithm', 'active-set');
			opts2 = optimset(opts2, 'MaxFunEvals', 50000, 'MaxIter', 50000, 'TolFun', 1e-012, 'TolX', 1e-012, 'FunValCheck','on');
			opts2 = optimset(opts2, 'LargeScale', 'off', 'Display', 'off', 'Algorithm', 'active-set');		
			try
				if ~exist('weight')
					[x1,fval1,eFlag1,output1] = fminsearch(@(param) sumLsqDirTuning(param, ...
								directions, y), x0, opts1);
				elseif exist('weight')
					[x1,fval1,eFlag1,output1] = fminsearch(@(param) sumLsqDirTuning(param, ...
								directions, y, weight), x0, opts1);	
				end
			catch
				if ~exist('weight')
					[x1,fval1,eFlag1,output1] = fmincon(@(param) sumLsqDirTuning(param, ...
								directions, y), x0, [], [], [], [], vlb, vub, [], opts2);
				elseif exist('weight')
					[x1,fval1,eFlag1,output1] = fmincon(@(param) sumLsqDirTuning(param, ...
								directions, y, weight), x0, [], [], [], [], vlb, vub, [], opts2);
				end
			end
			y_hat1 = circularNormalFunc(x1, directions);
			corrOne = corr(y, y_hat1);
			evOne = corrOne^2;

			ss0 = sumLsqDirTuning(x1, directions, y);
			ss0 = ss0/(length(y)-1);
			var0 = var(y);
			evVal = 1-ss0/var0;

			x = x1;
			fval = fval1;
			eFlag = eFlag1;
			output = evVal;
		catch
			x = [nan nan nan];
			fval = nan;
			eFlag = nan;
			output = nan;
		end
	end


function [x,fval,eFlag,output] = fitSumDirGaussian(initGuess, directions, y, weight)

	x0 = initGuess;
	
	vlb = [1 0 0.001 0];                	% lower bounds
	vub = [800 2*pi 2*pi 100];         		% upper bounds

	opts = optimset('fmincon');
	opts = optimset(opts, 'MaxFunEvals', 50000, 'MaxIter', 50000, 'TolFun', 1e-012, 'TolX', 1e-012, 'FunValCheck','on');
	opts = optimset(opts, 'LargeScale', 'off', 'Display', 'off', 'Algorithm', 'active-set');

	try
		if ~exist('weight')
			[x,fval,eFlag,output] = fmincon(@(param) sumLsqDirGaussian(param, directions, y), x0, [], [], [], [], vlb, vub, [], opts);
		elseif exist('weight')
			[x,fval,eFlag,output] = fmincon(@(param) sumLsqDirGaussian(param, directions, y, weight), x0, [], [], [], [], vlb, vub, [], opts);
		end
			
		y_hat = gaussianTunning(x, directions);

		ss0 = sumLsqDirGaussian(x, directions, y);
		ss0 = ss0/(length(y)-1);
		var0 = var(y);
		evVal = 1-ss0/var0;
		output = evVal;

	catch
		try
			opts1 = optimset('fminsearch');
			opts2 = optimset('fmincon');
			opts1 = optimset(opts1, 'MaxFunEvals', 50000, 'MaxIter', 50000, 'TolFun', 1e-012, 'TolX', 1e-012, 'FunValCheck','on');
			opts1 = optimset(opts1, 'LargeScale', 'off', 'Display', 'off', 'Algorithm', 'active-set');
			opts2 = optimset(opts2, 'MaxFunEvals', 50000, 'MaxIter', 50000, 'TolFun', 1e-012, 'TolX', 1e-012, 'FunValCheck','on');
			opts2 = optimset(opts2, 'LargeScale', 'off', 'Display', 'off', 'Algorithm', 'active-set');		
			try
				if ~exist('weight')
					[x1,fval1,eFlag1,output1] = fminsearch(@(param) sumLsqDirGaussian(param, ...
								directions, y), x0, opts1);
				elseif exist('weight')
					[x1,fval1,eFlag1,output1] = fminsearch(@(param) sumLsqDirGaussian(param, ...
								directions, y, weight), x0, opts1);	
				end
			catch
				if ~exist('weight')
					[x1,fval1,eFlag1,output1] = fmincon(@(param) sumLsqDirGaussian(param, ...
								directions, y), x0, [], [], [], [], vlb, vub, [], opts2);
				elseif exist('weight')
					[x1,fval1,eFlag1,output1] = fmincon(@(param) sumLsqDirGaussian(param, ...
								directions, y, weight), x0, [], [], [], [], vlb, vub, [], opts2);
				end
			end
			y_hat1 = gaussianTunning(x1, directions);
			corrOne = corr(y, y_hat1);
			evOne = corrOne^2;

			ss0 = sumLsqDirGaussian(x1, directions, y);
			ss0 = ss0/(length(y)-1);
			var0 = var(y);
			evVal = 1-ss0/var0;

			x = x1;
			fval = fval1;
			eFlag = eFlag1;
			output = evVal;
		catch
			x = [nan nan nan nan];
			fval = nan;
			eFlag = nan;
			output = nan;
		end
	end
	
function minValue = sumLsqDirTuning(param, directions, y, weight);

	y_hat = circularNormalFunc(param, directions);
	if nargin <= 3
		lsq_sum = (y - y_hat).^2;
	elseif nargin > 3
		lsq_sum = ((y - y_hat)./weight).^2;
	end
	
	minValue = sum(lsq_sum);

function minValue = sumLsqDirGaussian(param, directions, y, weight);

	y_hat = gaussianTunning(param, directions);
		
	if nargin <= 3
		lsq_sum = (y - y_hat).^2;
	elseif nargin > 3
		lsq_sum = ((y - y_hat)./weight).^2;
	end
	
	minValue = sum(lsq_sum);
	
function y_hat = gaussianTunning(param, directions);

	y_hat = param(1)*exp(-((directions - param(2)).^2)./(2*param(3)^2)) + param(4);
	
function y_hat = circularNormalFunc(param, directions);

	y_hat = param(1)*exp(param(2)*(cos(directions - param(3))-1));
	
	