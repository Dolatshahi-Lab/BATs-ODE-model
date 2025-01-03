function [passedtotalModelRuns, passedtotalParams, failedtotalModelRuns, failedtotalParams] = RB_eval_model_performance(calibrationInput, modelRuns, passCriteria)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%
% Evaluates model performance by comparing model run outcomes against
% experimental data provided by the user.
%
% INPUTS
%   - calibrationInput: structure that holds user input to calibration
%   program.
%   - modelRuns: structure that holds modelRuns output from the lhs_ode_run
%   function as well as other information such as paramMatrix.
%
% OUTPUTS
%   - passedtotalModelRuns -- a matrix that holds the runs that passedclear
%   following comparison to user-provided data.
%   - passedtotalParams:  a matrix holding the parameter values for the
%   runs that passed following comparison to user-provided data.
%   - failedtotalModelRuns:  a matrix that holds the runs that failed
%   following comparison to user-provided data. 
%   - failedtotalParams: a matrix that holds the parameter values for the
%   runs that failed following comparison to user-provided data. 
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

%name some tempvars for model Output:
TumorCancerOutput = modelRuns.odeOut{6,1}; % timesteps x number of runs for output 1 -- this is multipled by 1000 because the experimental data is in the thousands
TumorDonorOutput = modelRuns.odeOut{6,1}; % timesteps x number of runs for output 1 -- this is multipled by 1000 because the experimental data is in the thousands

% TumorCancerOutput = (1-(modelRuns.odeOut{1,1}./modelRuns.odeOut{11,1}))*100;
% TumorDonorOutput = (1-(modelRuns.odeOut{6,1}./modelRuns.odeOut{12,1}))*100;
% 

%name some tempvars for experimental data:
TumorCancerDS = calibrationInput.Tumor_cancerDS;
TumorDonorDS = calibrationInput.Tumor_donorDS;

% how does the output perform against dataset?
% [passedTumorRuns, runIndxTumorPassed, runIndxTumorFailed] = find_passed(TumorOutput, TumorDS,passCriteria);
% [passedBATsRuns, runIndDBATsPassed, runIndDBATsFailed] = find_passed(BATsOutput, BATsDS,passCriteria);
[~, runIndxTumorCancerPassed, runIndxTumorCancerFailed] = find_passed(TumorCancerOutput(modelRuns.odeSettings.time_index,:), TumorCancerDS,passCriteria);
[~, runIndxTumorDonorPassed, runIndxTumorDonorFailed] = find_passed(TumorDonorOutput(modelRuns.odeSettings.time_index,:), TumorDonorDS,passCriteria);

%find the columns that are in common between the two conditions (Tumor data
%and BATs data)
CommonIndxPassed = intersect(runIndxTumorCancerPassed,runIndxTumorDonorPassed);%, runIndxNBATsPassed);
% 
CommonIndxFailed = union(runIndxTumorCancerFailed, runIndxTumorDonorFailed);

%subset runs, LHS matrix, and param matrix to include only the ones that
%passed:
passedtotalParams = modelRuns.paramMatrix(CommonIndxPassed,:);
passedtotalModelRuns(:,:,1) = TumorCancerOutput(:,CommonIndxPassed);
passedtotalModelRuns(:,:,2) = TumorDonorOutput(:,CommonIndxPassed);

%subset runs and param matrix to also include those that did not pass:
failedtotalParams = modelRuns.paramMatrix(CommonIndxFailed,:);
failedtotalModelRuns(:,:,1) = TumorCancerOutput(:,CommonIndxFailed);
failedtotalModelRuns(:,:,2) = TumorDonorOutput(:,CommonIndxFailed);

assignin('base', 'failedtotalParams', failedtotalParams); %Not necessary long-term, helpful for debug right now tho
assignin('base', 'failedtotalModelRuns', failedtotalModelRuns); %Not necessary long-term, helpful for debug right now tho
assignin('base', 'passedtotalModelRuns_1', passedtotalModelRuns); %Not necessary long-term, helpful for debug right now tho
assignin('base', 'passedtotalParams_1', passedtotalParams); %Not necessary long-term, helpful for debug right now tho

end


function [passedModelRuns, RunIndxThatPassed, RunIndxThatFailed] = find_passed(modelOutput, expData, passCriteria)


%logical indexing to quickly check if the runs pass through the window of
%our data (time or divided by our passing criteria). 
LogicalArrayCondition1 = modelOutput(1:8,:) >= expData(1:8,1) / passCriteria;% [1;passCriteria;1;1;1;passCriteria;passCriteria;passCriteria;passCriteria;passCriteria];
LogicalArrayCondition2 = modelOutput(1:8,:) <= expData(1:8,2) * passCriteria;% [1;passCriteria;1;1;1;passCriteria;passCriteria;passCriteria;passCriteria;passCriteria];

LogicalArraytot = LogicalArrayCondition1 & LogicalArrayCondition2;

% cols with all ones --  meaning they passed the criteria
RunIndxThatPassed = find(all(LogicalArraytot==1)); 

%cols with all zeroes -- meaning they failed the criteria:
RunIndxThatFailed = find(~all(LogicalArraytot));

%use RunIndxThatPassed to pull those columns from the total runs matrix:
passedModelRuns = modelOutput(:,RunIndxThatPassed);
assignin('base', 'modelOutput', modelOutput); %Not necessary long-term, helpful for debug right now tho
assignin('base', 'expData', expData); %Not necessary long-term, helpful for debug right now tho
assignin('base', 'LogicalArraytot', LogicalArraytot); %Not necessary long-term, helpful for debug right now tho


end
