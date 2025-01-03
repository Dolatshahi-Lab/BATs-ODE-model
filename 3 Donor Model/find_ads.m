function [hdi, xdensityPass, xPass, xdensityFail, xFail, xdensitySub]=find_ads(passedtotalParamsVec, failedtotalParamsVec)
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% for a given set of parameter values, this file will find the density
% distribution of those parameter values, and calculate the ADS according
% to the algorithm set forth by Louis Joslyn in 2020. 
%
% Input: 
%    passedtotalParamsVec -- number of passed runs by 1 (the parameter that
%    is being calculated to find the HDR). 
%
%    failedtotalParamsVec -- number of failed runs by 1 (the parameter that
%    is being calculated to find the HDR, the same as the
%    passedtotalParamsVec parameter). 
%
%    coverage -- value that determines the region that we are searching
%    for. Value from 0-1, equal to value calibrationInput.portionParams.
%  
% Output:
%    hdi -- 1 by 2. The min and max parameter values that corespond to the
%    ADS identified region of parameter space where more pass runs exist
%    than fail runs.  These values will be used as the min and max for
%    sampling the parameter range on the next iteration of CaliPro.
%
%    xdensityPass -- 1 by 100. 1 by 100 set of density values that describe
%    the underlying pass dataset distribution.
%
%    xPass  -- the xvalues (parameter values) that correspond to the
%    xdensityPass values.  1 by 100.
%
%    xdensityFail -- 1 by 100. 1 by 100 set of density values that describe
%    the underlying fail dataset distribution.
%
%    xFail -- the xvalues (parameter values) that correspond to the
%    xdensityFail values.  1 by 100.
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


%find density estimation of the passed params
[xdensityPass,xPass] = ks_density_twice(passedtotalParamsVec);

%--------------------------------------------------------------------------
%find density of the failed params
%if there are no failedtotalParamsVec
if (isempty(failedtotalParamsVec) == 1 )
    xdensityFail = []; %make the density empty
    xFail = []; %make the xvalues that go with the density empty too 
else
    [xdensityFail,xFail] = ks_density_twice(failedtotalParamsVec);
end
%--------------------------------------------------------------------------
%only keep values from density that are included in the bounds of original
%parameter ranges:
indicesPass = find(xPass < min(passedtotalParamsVec) | xPass > max(passedtotalParamsVec));
xdensityPass(indicesPass) = 0;
%xPass(indicesPass) = [];

indicesFail = find(xFail < min(failedtotalParamsVec) | xFail > max(failedtotalParamsVec));
xdensityFail(indicesFail) = 0;
%xFail(indicesFail) = [];

%--------------------------------------------------------------------------
%subtract the failed from the passed:
xdensitySub = xdensityPass - xdensityFail;
%find all points along range where passed runs are greater than failed
%runs, so that we know the ranges where the runs (on the majority)
%performed well:
hdiRange = xPass(xdensitySub > 0);

%put the min and max into the highest density interval output.
hdi = [min(hdiRange), max(hdiRange)];

end