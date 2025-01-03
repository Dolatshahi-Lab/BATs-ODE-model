function RB_Calibration_Driver_pp(lhsODESettingsFileName, seedGenAndNum)
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%
% this is a Driver file for an automated, model-agnostic calibration
% strategy of complex systems that require unconventional objective
% function evaluations (i.e., non-continuous objective functions) and
% utilitzes an LHS sampling scheme for sweeping global parameter space.
% However, this process could be modified for any sampling scheme, and is
% not dependent on LHS sampling.  For example, Sobol sequencing could be
% used.
%
% CaliPro (as implemented below) is currently used with an ODE model that
% utilizes LHS sampling.  The sampling and model specification setup
% follows the process our lab has used for several years now, but the
% general approach of CaliPro should work for any model, including
% stochastic models like ABMs, and we hope to generalize this framework for
% future use.
%
% The steps of the calibration method are the following:
%       1. Provide the Program with the Necessary Information to Run.
%            - initially, this includes the model itself, the starting
%              parameter ranges, any data to compare against (experimental
%              data, for example) and the termination condition.
%       2. Perform Sampling According to the Sampling Scheme (LHS). 
%       3. Run the Model 
%              - steps 2 and 3 might be combined in using the
%              Predator-Prey model system evaluated using our lhs_ode lab
%              setup
%        4. Evaluate the Performance of the Model Against Experimental Data.
%            - Create Passed Model Runs Data structure
%            - Create Passed Model Runs Params Data structure
%            - if performace satisfies termination condition, then leave 
%              the calibration program.
%       5. Perform Density (Distribution) Comparisons between Passed Model
%          Parameter and Total Model Parameter Ranges.
%            - Evaluate the densest portions
%              of the passed parameter set ranges. 
%            - The range of these densest regions are the new parameter
%              ranges for the next iteration of the calibration program. 
%            - AlternateRun designates if the user wishes to use the
%            alternate method of comparing passed runs vs failed runs via
%            substitution of the latter from the former.  This substitution
%            gives a range of parameter values wherein the runs mostly
%            passed. 
% 
% Repeat this process (particularly steps 2-5) until condition is
% satisfied. As an arbitrary example, the condition might be satisfied when
% 90% of individiuals have an infection where all granulomas can control
% the bacterial growth. 
% 
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
%     1. Provide the Program with the Necessary Information to Run.
% -------------------------------------------------------------------------

% model information is stored in the model settings file 
% Load the settings file - everything returned in the odeSettings structure.
odeSettings = lhs_ode_get_run_settings_new(lhsODESettingsFileName);

% if seed information was input, then assign that information to the seed
% generator and seed number.
if exist('seedGenAndNum', 'var')
    rng(seedGenAndNum);
end
% save the current generator settings for reproduciblity.
% this approach allows for long-term reproducibility by saving not only the
% seed, but also the seed generator.
odeSettings.seedGenSettings = rng;

% eventually, load the calibration settings file - everything returned in
% the calibration structure. For now, just put the information inside the
% driver file

calibrationInput.HDR = false;

%How long should this run for unsuccesfully?
calibrationInput.maxIter = 5;

%What is the percentage of runs that need to satisfy our criteria for the
%termination condition to be executed?
calibrationInput.terminate = 0.7;

%What is the portion of parameters that we should be selecting
%within the parameter ranges at the end of each iteration? 
calibrationInput.portionParams = 0.5; 

%how flexible is the passing criteria? 
calibrationInput.passCriteria = 2;

%what are the names of the experimental datasets?  Datasets to compare
%against model outcomes  
calibrationInput.MCF7DS = csvread('TMinMax_MCF7.csv');
calibrationInput.CAMA1DS = csvread('TMinMax_CAMA1.csv');
calibrationInput.T47DDS = csvread('TMinMax_T47D.csv');


calibrationInput.passCriteriaInitial = 2.5;
calibrationInput.passCriteriaNonInitial = 1.75;

%set up the array that holds all of the calibration results for each run:
CalibrationOutputArray = {};
%loop through, until condition is satisfied:
runIteration = 1;
terminateCalib = false; 
while ~terminateCalib
    %for runIteration = 1:1
    % -------------------------------------------------------------------------
    %     2. Perform Sampling According to the Sampling Scheme (LHS). 
    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------
    %     3. Run the Model.
    % -------------------------------------------------------------------------
    %for this model framework, LHS sampling and model execution happen through
    %calling the lhs_ode_run_new_c function.  the lhs_ode_run_new function has
    %been altered to work better within the calibration framework.
    modelRuns = RB_lhs_ode_run_new_c(odeSettings,runIteration);


    % -------------------------------------------------------------------------
    %     4. Evaluate the Performance of the Model Against Experimental Data.
    % -------------------------------------------------------------------------

    %make fit criteria stricter after initial run:
    passCriteria = calibrationInput.passCriteriaInitial;
    if runIteration > 2
        passCriteria = calibrationInput.passCriteriaNonInitial;
    end

    % Compare model output against the experimental data
    [passedtotalModelRuns, passedtotalParams, failedtotalModelRuns, failedtotalParams] = RB_eval_model_performance(calibrationInput, modelRuns, passCriteria);
    % save output/plot results
    modelRuns.passedtotalModelRuns = passedtotalModelRuns;
    modelRuns.passedtotalParams = passedtotalParams;
    modelRuns.failedtotalModelRuns = failedtotalModelRuns;
    modelRuns.failedtotalParams = failedtotalParams;

    % evaluate the performance according to termination conditions
    %leave loop once condition is satisfied:
    if (size(passedtotalModelRuns,2) / odeSettings.NR) > calibrationInput.terminate
        disp('Calibration Complete According to User Specifications')
        CalibrationOutputArray{runIteration} = modelRuns;
        terminateCalib = true;
        break
    end 

    % -------------------------------------------------------------------------
    %     5. Perform Density (Distribution) Comparisons between Passed Model
    %        Parameter and Total Model Parameter Ranges.
    % -------------------------------------------------------------------------   
    %   if user specfied,             
    if ~(calibrationInput.HDR) %find the param ranges using alternate method via subtraction method
        [hdr_param1, density1, xi1, xdensityFail1, xFail1, xdensitySub1] = find_ads(passedtotalParams(:,1), failedtotalParams(:,1));
        [hdr_param2, density2, xi2, xdensityFail2, xFail2, xdensitySub2] = find_ads(passedtotalParams(:,2), failedtotalParams(:,2));
        [hdr_param3, density3, xi3, xdensityFail3, xFail3, xdensitySub3] = find_ads(passedtotalParams(:,3), failedtotalParams(:,3));
        [hdr_param4, density4, xi4, xdensityFail4, xFail4, xdensitySub4] = find_ads(passedtotalParams(:,4), failedtotalParams(:,4));
        % [hdr_param5, density5, xi5, xdensityFail5, xFail5, xdensitySub5] = find_ads(passedtotalParams(:,5), failedtotalParams(:,5));
        % [hdr_param6, density6, xi6, xdensityFail6, xFail6, xdensitySub6] = find_ads(passedtotalParams(:,6), failedtotalParams(:,6));
               
    else %find hdr for each param in passedtotal params. Return failed info too   
        [hdr_param1, density1, xi1, xdensityFail1, xFail1] = find_hdr(passedtotalParams(:,1), failedtotalParams(:,1), calibrationInput.portionParams);
        [hdr_param2, density2, xi2, xdensityFail2, xFail2] = find_hdr(passedtotalParams(:,2), failedtotalParams(:,2), calibrationInput.portionParams);
        [hdr_param3, density3, xi3, xdensityFail3, xFail3] = find_hdr(passedtotalParams(:,3), failedtotalParams(:,3), calibrationInput.portionParams);
        [hdr_param4, density4, xi4, xdensityFail4, xFail4] = find_hdr(passedtotalParams(:,4), failedtotalParams(:,4), calibrationInput.portionParams);
        % [hdr_param5, density5, xi5, xdensityFail5, xFail5] = find_hdr(passedtotalParams(:,5), failedtotalParams(:,5), calibrationInput.portionParams);
        % [hdr_param6, density6, xi6, xdensityFail6, xFail6] = find_hdr(passedtotalParams(:,6), failedtotalParams(:,6), calibrationInput.portionParams);
          
    end

    %load the new parameter ranges, based on the output of the find_hdr or
    %find_hdr_alt.  If the number of parameters is larger, there could be some
    %sort of function or (at the very least) a for loop that does this for all
    %the parameters. Could also include a for loop for the if else statement
    %above.
    odeSettings.parameters = ...
        {
        { 'aKill_N', 'lu', abs(hdr_param1(1)), abs(hdr_param1(2))} ...
        { 'aKill_P', 'lu', abs(hdr_param2(1)), abs(hdr_param2(2))} ...
        { 'aKill_L', 'lu', abs(hdr_param3(1)), abs(hdr_param3(2))} ...
        { 'aKill_D', 'lu', abs(hdr_param4(1)), abs(hdr_param4(2))} ...
        
        };

    %update modelRuns structure with the Density values for those parameters
    %whose simulation resulted in a pass/fail, the parameter range for those
    %runs and, if specified, the density values that represent the subtraction
    %of passed parameters from failed parameters
    modelRuns.ParamsDensityPass = [density1; density2; density3; density4];%; density5; density6];%; density7; density8; density9; density10; density11; density12; density13; density14; density15; density16; density17; density18; density19; density20; density21; density22; density23];
    modelRuns.ParamsDensityRangePass = [xi1; xi2; xi3; xi4];%; xi5; xi6];%; xi7; xi8; xi9; xi10; xi11; xi12; xi13; xi14; xi15; xi16; xi17; xi18; xi19; xi20; xi21; xi22; xi23];
    modelRuns.ParamsDensityFail = [xdensityFail1; xdensityFail2; xdensityFail3; xdensityFail4];%; xdensityFail5; xdensityFail6];%; xdensityFail7; xdensityFail8; xdensityFail9; xdensityFail10; xdensityFail11; xdensityFail12; xdensityFail13; xdensityFail14; xdensityFail15; xdensityFail16; xdensityFail17; xdensityFail18; xdensityFail19; xdensityFail20; xdensityFail21; xdensityFail22; xdensityFail23];
    modelRuns.ParamsDensityRangeFail = [xFail1; xFail2; xFail3; xFail4];%; xFail5; xFail6];%; xFail7; xFail8; xFail9; xFail10; xFail11; xFail12; xFail13; xFail14; xFail15; xFail16; xFail17; xFail18; xFail19; xFail20; xFail21; xFail22; xFail23];
    if ~(calibrationInput.HDR) 
        modelRuns.DensitySub = [xdensitySub1; xdensitySub2; xdensitySub3; xdensitySub4];%; xdensitySub5; xdensitySub6];%; xdensitySub7; xdensitySub8; xdensitySub9; xdensitySub10; xdensitySub11; xdensitySub12; xdensitySub13; xdensitySub14; xdensitySub15; xdensitySub16; xdensitySub17; xdensitySub18; xdensitySub19; xdensitySub20; xdensitySub21; xdensitySub22; xdensitySub23];
    end
    modelRuns.hdrRange = [hdr_param1; hdr_param2; hdr_param3; hdr_param4];%; hdr_param5; hdr_param6];%; hdr_param7; hdr_param8; hdr_param9; hdr_param10; hdr_param11; hdr_param12; hdr_param13; hdr_param14; hdr_param15; hdr_param16; hdr_param17; hdr_param18; hdr_param19; hdr_param20; hdr_param21; hdr_param22; hdr_param23];

    %update CalibrationOutputArray with all the information from this run
    %iteration in the calibration program. 
    CalibrationOutputArray{runIteration} = modelRuns;

    %update the runIteration counter
    runIteration = runIteration + 1;

    % evaluate the performance according to termination conditions
    %leave loop once condition is satisfied:
    if (runIteration > calibrationInput.maxIter)
        disp('Calibration Met Max Iteration Count. Evaluate initial parameter ranges or model construction')
        terminateCalib = 1;
    end 

end  %end while loop that performs the calibration

% Assign variables to the base workspace so they will be saved after
% leaving the calibration.
assignin('base', 'CalibrationOutputArray', CalibrationOutputArray);
assignin('base', 'calibrationInput', calibrationInput);
assignin('base', 'passedtotalModelRuns', passedtotalModelRuns); %Not necessary long-term, helpful for debug right now tho
assignin('base', 'failedtotalModelRuns', failedtotalModelRuns); %Not necessary long-term, helpful for debug right now tho
assignin('base', 'passedtotalParams', passedtotalParams); %Not necessary long-term, helpful for debug right now tho
assignin('base', 'odeSettings', odeSettings); %Not necessary long-term, helpful for debug right now tho
end %end calibration function
