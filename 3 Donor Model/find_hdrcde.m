function [hdr, den, x, xdensityFail, xFail, falpha] = find_hdrcde(passedtotalParamsVec, failedtotalParamsVec, prob)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% This function is a Matlab implementation of Rob Hyndman's Highest Density
% Region (HDR) R package titled "hdrcde".  It does not include the full
% functionality of that package.  Namely, it does not allow for multiple
% probabilities entered as an input (you can only provide one probability)
% and it does not output the distribution's mode, density plot, or boxplot.
% 
% The algorithm is outlined in Hyndman's 1996 article in The American
% Statistician called "Computing and Graphing Highest Density Regions". In
% it, he argues that HDR's are the most appropriate method for summarizing
% a probability distribution. The HDR can find the smallest interval of x
% in which the probability distribution, f(x), is densest, including the
% possibility of a disjoint interval, should the distribution be
% multimodal. 
%
% The user can read the paper for more information, but briefly, HDR works
% in the following manner. Let f(x) be the density function of a random
% variable X. Then the 100(1 - a)% HDR is the subset R(falpha) of the
% sample space of X such that R(falpha) = {x: f(x) >= falpha} where falpha
% is the largest constant such that Pr(X is an element of 
% R(falpha)) >= 1 - alpha. To be clear, this approach avoids explicit
% integration under f(x) by providing a density quantile approach.
%
% 
%
% Author: Louis Joslyn 
% Date: May 10, 2019
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------


%find the density estimation using the build-in matlab estimator
[xdensityPass,xPass, bandwidth1_P] = ksdensity(passedtotalParamsVec);
% run density estimate again, but with the bandwidth pre-specified to be
% 1/2 of what the embedded selection process of the bandwidth is within the
% Matlab ksdensity function.  Empirical research shows that, for whatever
% reason, the bandwidth selected by this algorithm is too smooth.  Dividing
% that number by 1/2 makes the density approximation a much more realistic
% representation of the underlying data distributon.
[den,x] = ksdensity(passedtotalParamsVec, 'bandwidth', bandwidth1_P/2);

% identify alpha1.  alpha1 = 1-alpha. where alpha is the user inputed
% probability. 
alpha1 = 1-(prob/100);
% calculate falpha
[falpha, fx] = calc_falpha(x,den,passedtotalParamsVec,alpha1);
%find HDR using falpha
[hdr,falpha] = hdr_ends(x, den,falpha);



%find the density for the failed parameters for plotting purposes following
% calibration completition.

if (isempty(failedtotalParamsVec) == 1 ) %if there are no failedtotalParamsVec
    xdensityFail = []; %make the density empty
    xFail = []; %make the xvalues that go with the density empty too
    
else %find the densities for the failed parameters as well    
    %find the density estimation using the build-in matlab estimator
    [xdensityFail,xFail, bandwidth1_F] = ksdensity(failedtotalParamsVec);
    %run again, but with the bandwidth prespecified to be 1/2 of what the
    %embedded selection process of the bandwidth is
    % for the unedited version.  Empirical research done with Marissa Renardy
    % shows that, for whatever reason, the bandwidth selected by this algorithm
    % is too smooth.  Dividing that number by 1/2 makes the density
    % approximation a much more realistic representation of the underlying data
    % distributon.
    [xdensityFail,xFail] = ksdensity(failedtotalParamsVec, 'bandwidth', bandwidth1_F/2);
end




end
    

function [hdr, falpha] = hdr_ends(x, den, falpha)
%find length of x
n = length(x);
%if falpha is larger than your density, you will not be able to find hdr of
%the distribution, so return nothing. 
if (falpha > max(den))
    hdr = NaN;
else %okay, so density is larger than falpha, good, you can proceed.
    dd = den - falpha; %subtract falpha from the density to get a sense of what is larger than falpha
    ddUP = dd(2:n).*dd(1:(n-1)); %create array to represent dd ranges
    indexAll = (1:n-1); %create array of same size as ddUP
    %subset array by only those values less than or equal to 0. these
    %should be the points where, along the range of den, the den-falpha
    %becomes negative. 
    indexSub = indexAll(ddUP<=0); 
    ni = length(indexSub); 
    %create intercept, a matrix of size length(indexSub).
    intercept = zeros(1, ni);
    
    if (ni>0)
        for j=1:ni
            idx = [indexSub(j),indexSub(j)+1];
            intercept(j) = interp1(den(idx), x(idx), falpha);
        end 
    
    end 
    intercept = unique(intercept); %unique function sorts returns sorted list of unique values
    ni = length(intercept);
    if (ni == 0)
        intercept = [x(1),x(n)];
    end
    x1 = 0.5*(intercept(1) + x(1));
    x2 = 0.5*(intercept(ni) + x(n));
    if(interp1(x,den,x1) > falpha)
      intercept = [0,intercept];
    end
    if(interp1(x,den,x2) > falpha)
        intercept = [intercept,0];
    end
    hdr=intercept;
end
end




function [falpha, fx] = calc_falpha(x, den, value, alpha1)
%      Calculates falpha needed to compute HDR of density den.
%      Input: den = density on grid.
%               x = independent observations on den
%           alpha = level of HDR
%      Output: falpha = The value of the density at the boundaries of each HDR.
%                  fx = interpolated density

%interpolate across the density to provide continuous points of evaluation
%and call it fx. 
fx = interp1(x,den,value, 'linear');
%falpha is the the alpha quantile of f(x), the probability density
%function. falpha can be estimated as a sample quantile from a set of iid
%random variables with the same distribution as f(x) --- thats what the
%input, fx, is. 
falpha = quantile(sort(fx), alpha1);

end 
