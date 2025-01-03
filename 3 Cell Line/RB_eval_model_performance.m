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
MCF7Output = modelRuns.odeOut{1,1}; % timesteps x number of runs for output 1 -- this is multipled by 1000 because the experimental data is in the thousands
CAMA1Output = modelRuns.odeOut{6,1}; % timesteps x number of runs for output 1 -- this is multipled by 1000 because the experimental data is in the thousands
T47DOutput = modelRuns.odeOut{11,1}; % timesteps x number of runs for output 1 -- this is multipled by 1000 because the experimental data is in the thousands

%name some tempvars for experimental data:
MCF7DS = calibrationInput.MCF7DS;
CAMA1DS = calibrationInput.CAMA1DS;
T47DDS = calibrationInput.T47DDS;

% how does the output perform against dataset?
[passedMCF7Runs, runIndxMCF7Passed, runIndxMCF7Failed] = find_passed(MCF7Output(modelRuns.odeSettings.time_index,:), MCF7DS,passCriteria);
[passedCAMA1Runs, runIndxCAMA1Passed, runIndxCAMA1Failed] = find_passed(CAMA1Output(modelRuns.odeSettings.time_index,:), CAMA1DS,passCriteria);
[passedT47DRuns, runIndxT47DPassed, runIndxT47DFailed] = find_passed(T47DOutput(modelRuns.odeSettings.time_index,:), T47DDS,passCriteria);


%find the columns that are in common between the two conditions (Tumor data
%and BATs data)
CommonIndxPassed = intersect(runIndxMCF7Passed,runIndxCAMA1Passed);
CommonIndxPassed = intersect(CommonIndxPassed, runIndxT47DPassed);
% 
CommonIndxFailed = union(runIndxMCF7Failed, runIndxCAMA1Failed);
CommonIndxFailed = union(CommonIndxFailed, runIndxT47DFailed);

%subset runs, LHS matrix, and param matrix to include only the ones that
%passed:
passedtotalParams = modelRuns.paramMatrix(CommonIndxPassed,:);
passedtotalModelRuns(:,:,1) = MCF7Output(:,CommonIndxPassed);
passedtotalModelRuns(:,:,2) = CAMA1Output(:,CommonIndxPassed);
passedtotalModelRuns(:,:,3) = T47DOutput(:,CommonIndxPassed);

%subset runs and param matrix to also include those that did not pass:
failedtotalParams = modelRuns.paramMatrix(CommonIndxFailed,:);
failedtotalModelRuns(:,:,1) = MCF7Output(:,CommonIndxFailed);
failedtotalModelRuns(:,:,2) = CAMA1Output(:,CommonIndxFailed);
failedtotalModelRuns(:,:,3) = T47DOutput(:,CommonIndxFailed);

assignin('base', 'failedtotalParams', failedtotalParams); %Not necessary long-term, helpful for debug right now tho
assignin('base', 'failedtotalModelRuns', failedtotalModelRuns); %Not necessary long-term, helpful for debug right now tho
assignin('base', 'passedtotalModelRuns_1', passedtotalModelRuns); %Not necessary long-term, helpful for debug right now tho
assignin('base', 'passedtotalParams_1', passedtotalParams); %Not necessary long-term, helpful for debug right now tho
assignin('base', 'passedtotalParams', passedtotalParams); %Not necessary long-term, helpful for debug right now tho
assignin('base', 'modelRuns_odeOut', modelRuns.odeOut); %Not necessary long-term, helpful for debug right now tho

end


function [passedModelRuns, RunIndxThatPassed, RunIndxThatFailed] = find_passed(modelOutput, expData, passCriteria)


%logical indexing to quickly check if the runs pass through the window of
%our data (time or divided by our passing criteria). 
LogicalArrayCondition1 = modelOutput([1,3:8],:) >= expData([1,3:8],1) / passCriteria;
LogicalArrayCondition2 = modelOutput([1,3:8],:) <= expData([1,3:8],2) * passCriteria;
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
