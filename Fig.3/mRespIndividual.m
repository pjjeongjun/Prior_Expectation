function R = mRespIndividual(gain, direction, PD, sigD, sponR, cparams1, cparams2, cparams3, cparams4)
% gaussian shape tuning curve
% direction = stimulus dir
% sigD = tuning width

% pd range: -180 ~ 180
if PD > 180
	PD = PD - 360;
end
if PD < -180
	PD = PD + 360;
end

% degree to radian
dirD = (direction - PD)*pi/180.0;
sigD = sigD*pi/180.0; 

% Gaussian
if isnan(cparams1(1)) && isnan(cparams2(1)) && isnan(cparams3(1)) && isnan(cparams4(1))
    R = gain*exp(-dirD.^2./(2*sigD^2))+sponR;
elseif ~isnan(cparams1(1)) && ~isnan(cparams2(1)) && ~isnan(cparams3(1)) && ~isnan(cparams4(1))
    R = cparams1*exp(-((direction*pi/180.0 - cparams2).^2)./(2*cparams3^2)) + cparams4;
%     R = cparams1*exp(-(dirD.^2)./(2*cparams3^2)) + cparams4;
    
% Circular Gaussian
elseif ~isnan(cparams1(1)) && ~isnan(cparams2(1)) && ~isnan(cparams3(1)) && isnan(cparams4(1))
%     R = cparams1*exp(cparams2*(cos(direction*pi/180.0 - cparams3)-1));
    R = cparams1*exp(cparams2*(cos(dirD)-1)); % Circular Gaussian
end
  