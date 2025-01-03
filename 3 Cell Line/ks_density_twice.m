function [xDen, xVal] = ks_density_twice(paramVec)
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% ks_density_twice.m is a function that calls the ksdensity function twice
% with a slight change in bandwidth the second time.  We found, through
% empirical research of comparing distribution histograms against the
% built-in matlab ks density estimator, that  we needed to run the function
% twice, but with the bandwidth prespecified to be 1/2 of what the embedded
% selection process of the bandwidth is for the unedited version. 
% 
% input:
%    paramVec -- a vector of parameter values that correspond to either the
%    pass or fail set of runs. (number of varied parameters by number of
%    runs that satisfied the criteria)
% output:
%    xDen -- 1 by 100 set of density values that describe the underlying
%    data distribution.
%
%    xVal -- the xvalues (parameter values) that correspond to the xDen
%    values.  1 by 100. 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%use built-in matlab function, ksdensity to create a kernel based density
%estimation of the underlying data distribution.
[xDen1,xVal1, bandwidth1_P] = ksdensity(paramVec);
%run again, but with the bandwidth prespecified to be 1/2 of what the
%embedded selection process of the bandwidth is
% for the unedited version.  Empirical research done with Dr. Marissa
% Renardy shows that, for whatever reason, the bandwidth selected by this
% algorithm is too smooth.  Dividing that number by 1/2 makes the density
% approximation a much more realistic representation of the underlying data
% distributon.
[xDen,xVal] = ksdensity(paramVec, 'bandwidth', bandwidth1_P/2);


end