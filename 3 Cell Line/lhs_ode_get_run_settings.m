function ode = lhs_ode_get_run_settings(lhsSettingsFileName)

% Get the ODE LHS settings from a settings file.
%
% The caller will perform checks on some of the settings, since
% different callers will have different settings requirements in some
% cases.

if isempty(lhsSettingsFileName)
    title = 'Select ODE LHS settings file';
    default = 'lhs_ode_settings.m';
    [lhsSettingsFileName, lhsDirectory, ~] = uigetfile('*.m', title, default);
    
    % The eval function (used below) expects a string, not a file name.
    % A file name only works with eval if it does not contain any path
    % information, does not have the trailing '.m' and is in the Matlab path.
    lhsSettingsPathName = [lhsDirectory, filesep, lhsSettingsFileName];

else
    % The directory is not a separate item if specified as a function
    % argument.
    lhsDirectory = '';
    lhsSettingsPathName = lhsSettingsFileName;
end    

if lhsSettingsFileName == 0
    msgIdent = 'LHS:NoSettingsFile';
    msg = 'No ODE LHS settings file specified.\n';
    throw(MException(msgIdent, msg));
end

% Load the LHS settings.
try
    if ~isempty(lhsDirectory)
        curDir = cd(lhsDirectory);
    end
    [~, name, ~] = fileparts(lhsSettingsFileName);
    tmpStr = strcat('ode = ',name, ';');
    eval(tmpStr);
    if ~isempty(lhsDirectory)
        cd(curDir);
    end
catch e
    msgIdent = 'LHS:ODE:SettingsFileEvalError';
    msg = ['Unable to load ODE LHS settings from file "%s".\n%s', ...
            '\nThis might be due to a typo or', ...
            ' stray characters in the settings file.', ...
            '\nThe settings file must contain valid Matlab code.'];
    throw(MException(msgIdent, msg, lhsSettingsPathName, e.message));
end

assignin('base','lhsOdeSettingsFileName', lhsSettingsFileName);
assignin('base', 'lhsDirectory', lhsDirectory);

% Perform some consistency checks on the loaded settings.

if (~exist('ode', 'var'))
    msgIdent = 'LHS:ODE:MissingSetting';
    msg = 'The ode parameter structure does not exist.';
    throw(MException(msgIdent, msg));    
end

if(isfield(ode, 'y0') && isempty(ode.y0))
    msg = 'The y0 (ODE initial conditions) setting does not exist or is empty.';
    msgIdent = 'LHS:ODE:MissingSetting';
    throw(MException(msgIdent, msg));
end

if (isempty(ode.tspan))
    msgIdent = 'LHS:ODE:MissingSetting';
    msg = 'The tspan (model output time points) setting does not exist.';
    throw(MException(msgIdent, msg));
end
% retrieve tspan from workspace
if (~isnumeric(ode.tspan))
    msgIdent = 'LHS:ODE:LHS:InvalidTspanList';
    msg = 'The tspan list of model time points is not numeric.';
    throw(MException(msgIdent, msg));  
end

% tspan must be a montonically increasing sequence of numeric values.
% The values in tspan will have a Matlab type of double.
t_previous = -1;
for t = ode.tspan
    if t < 0
        msgIdent = 'LHS:ODE:InvalidTimePoint';
        msgfmt = 'The tspan time point value of %d is < 1';
        msg = sprintf(msgfmt, t);
        throw(MException(msgIdent, msg));       
    elseif t <= t_previous
        msgIdent = 'LHS:ODE:InvalidTimePoint';
        msgfmt = ['The tspan time point value of %d' ...
                    ' is <= the previous time point of %d'];
        msg = sprintf(msgfmt, t, t_previous);
        throw(MException(msgIdent, msg));                   
    end
    
    t_previous = t;
end

% time_points must be a montonically increasing sequence of numeric values.
% The values in time_points will have a Matlab type of double.
% Each element of time points must also be an element of tspan.
% Also define the indices into span corrsponding to each time_point
% element. These will be used to index into result matrices which will have
% one row for each element in tspan.
%
% For example, suppose tspan = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0], and
% time_points = [1.5, 2.5] then time_index = [3, 5], since 3 is the index
% in tspan for 1.5, and 5 is the index in tspan for 2.5. Result matrices
% would have 6 rows, one for each element of tspan. If we want to perform
% analysis on the rows corresponding to time points 1.5 and 2.5, that
% corresponds to using rows 3 and 5 of those result matrices.
t_previous = -1;
ode.time_index = zeros(1, length(ode.time_points));
i = 1; % The next element of ode.time_index to be defined.
for t = ode.time_points
    if t < 0
        msgIdent = 'LHS:ODE:InvalidTimePoint';
        msgfmt = 'The time_points time point value of %d is < 1';
        msg = sprintf(msgfmt, t);
        throw(MException(msgIdent, msg));       
    elseif t <= t_previous
        msgIdent = 'LHS:ODE:InvalidTimePoint';
        msgfmt = ['The time_points time point value of %d' ...
                    ' is <= the previous time point of %d'];
        msg = sprintf(msgfmt, t, t_previous);
        throw(MException(msgIdent, msg));                   
    end
    
    [flag, tspan_index] = ismember(t, ode.tspan);
    if flag == 1
        ode.time_index(i) = tspan_index;
    else
        msgIdent = 'LHS:ODE:InvalidTimePoint';
        msgfmt = 'The time_points time point value of %d is not in tspan';
        msg = sprintf(msgfmt, t);
        throw(MException(msgIdent, msg));
    end
    t_previous = t;
    i = i + 1;
end

% Check the LHS settings needed for LHS setup.
varyingParameterCount = length(ode.PRCC_var);
if varyingParameterCount == 0
    % Not varying any parameters. Should only do 1 run in this case.
    is_numeric_min(ode.NR, 1, 'NR', 'LHS:ODE:RUN:InvalidSetting');
    if (ode.NR ~= 1)
        msgIdent = 'LHS:ODE:RUN:InvalidSetting';
        msgfmt = ['There are no varying parameters but the number of runs' ...
                  ' (NR) of %d is not equal to 1. It does not make sense' ...
                  ' to do more than one run if not varying any parameters.'];
        msg = sprintf(msgfmt, ode.NR);
        throw(MException(msgIdent, msg));                   
    end
elseif varyingParameterCount == 1
    % Not useful to vary only 1 parameter - treat as an error.
    msgIdent = 'LHS:ODE:PRCC_var:InvalidSetting';
    msg = ['There is only 1 varying parameter. PRCC cannot be performed' ...
              ' on only 1 varying parameter (i.e. for 2 or more runs) and' ...
              ' it does not make sense to vary 1 parameter for only 1 run.'];
    throw(MException(msgIdent, msg));
else
    % Varying more than 1 parameter, must do at least 2 runs.
    is_numeric_min(ode.NR, 2, 'NR', 'LHS:ODE:RUN:InvalidSetting');
end

if (~isa(ode.ode_model_handle, 'function_handle'))
    msgIdent = 'LHS:ODE:RUN:MissingSetting';
    msg = 'The ODE model function handle is not defined.';
    throw(MException(msgIdent, msg));
end

% Check that each element of PRCC_var is also a member of PRCC_labels.
badVarList = {};
for i = 1:length(ode.PRCC_var)
    p = ode.PRCC_var{i};
    if isempty(intersect(p, ode.PRCC_labels))
        badVarList{end+1} = p;
    end
end

if ~isempty(badVarList)
    fprintf('The following elements of PRCC_var are not in PRCC_labels.\n');
    fprintf('Check for mispelled names.\n');
    for i = 1:length(badVarList)
        fprintf('%s ', badVarList{i});
    end
    fprintf('\n');
    msgIdent = 'LHS:ODE:PRCC_var:NotExist';
    throw(MException(msgIdent, ''));
end

% Check the pulse items, if defined.
if isfield(ode, 'pulseTimesteps')

    % pulseTimesteps should not be empty.
    if isempty(ode.pulseTimesteps)
        msgIdent = 'LHS:ODE:InvalidPulseTimesteps';
        msgfmt = 'The pulseTimesteps vector is empty';
        msg = sprintf(msgfmt, t);
        throw(MException(msgIdent, msg));       
    end
    
    % pulseTimesteps must be a montonically increasing sequence of numeric
    % values. The values will have a Matlab type of double. 
    % Also each pulse time point should be in the tspan vector.
    pts_previous = -1;
    ode.pulseTimestepIdx = [];
    for pts = ode.pulseTimesteps
        if pts < 0
            msgIdent = 'LHS:ODE:InvalidPulseTimePoint';
            msgfmt = 'The pulse time point value of %d is < 0';
            msg = sprintf(msgfmt, pts);
            throw(MException(msgIdent, msg));       
        elseif pts <= pts_previous
            msgIdent = 'LHS:ODE:InvalidPulseTimePoint';
            msgfmt = ['The pulse time point value of %d' ...
                        ' is <= the previous time point of %d'];
            msg = sprintf(msgfmt, pts, pts_previous);
            throw(MException(msgIdent, msg));                   
        end

        pulseTimestepIdx = find(ode.tspan == pts);
        if isempty(pulseTimestepIdx)
            msgIdent = 'LHS:ODE:SETTINGS:BadPulseTimestep';
            msg = 'pulseTimestep of  %d not found in tspan\n';
            error(msgIdent, msg, pts);
        else
            if pulseTimestepIdx == 1
                msgIdent = 'LHS:ODE:SETTINGS:BadPulseTimestepIdx';
                msg = 'pulseTimestepIdx is 1. Can''t pulse on the 1st time step.\n';
                error(msgIdent, msg);
            end
        end
        
        ode.pulseTimestepIdx = [ode.pulseTimestepIdx pulseTimestepIdx];
        
        pts_previous = pts;
    end
    
    if ode.pulseIdx < 1
        msgIdent = 'LHS:ODE:SETTINGS:BadPulseIdx';
        msg = 'pulseIdx of %d is < 1. It must be in the range of initial conditions, 1 to %d.\n';
        error(msgIdent, msg, ode.pulseIdx, ode.initialConditionCount);
    end
       
    if ode.pulseIdx > ode.initialConditionCount
        msgIdent = 'LHS:ODE:SETTINGS:BadPulseIdx';
        msg = 'pulseIdx of %d is > the number of initial conditions. It must be in the range of initial conditions, 1 to %d.\n';
        error(msgIdent, msg, ode.pulseIdx, ode.initialConditionCount);
    end
    
    if isfield(ode, 'pulseAmountIdx')
        
        % pulseAmountIdx should not be empty.
        if isempty(ode.pulseAmountIdx)
            msgIdent = 'LHS:ODE:InvalidPulseAmountIdx';
            msgfmt = 'The pulseAmountIdx vector is empty';
            msg = sprintf(msgfmt, t);
            throw(MException(msgIdent, msg));       
        end

        if length(ode.pulseAmountIdx) ~= length(ode.pulseTimesteps)
            msgIdent = 'LHS:ODE:InvalidPulseAmountIdx';
            msgfmt = 'pulseAmountIdx, %d elements, is not the same length as ode.pulseTimesteps, %d elements.';
            msg = sprintf(msgfmt, length(ode.pulseAmountIdx), length(ode.pulseTimesteps));
            throw(MException(msgIdent, msg));       
        end
        
        
        idx_previous = -1;
        for idx = ode.pulseAmountIdx
            if idx < 0
                msgIdent = 'LHS:ODE:SETTINGS:BadPulseAmountIdx';
                msg = 'pulseAmountIdx of %d is < 1. It must be in the range of initial conditions, 1 to %d.\n';
                error(msgIdent, msg, idx, ode.initialConditionCount);
            elseif idx <= idx_previous
                msgIdent = 'LHS:ODE:InvalidPulseAmountIdx';
                msgfmt = ['The pulse amount index of %d' ...
                            ' is <= the previous one of %d'];
                msg = sprintf(msgfmt, idx, idx_previous);
                throw(MException(msgIdent, msg));
            elseif idx > ode.initialConditionCount
                msgIdent = 'LHS:ODE:SETTINGS:BadPulseAmountIdx';
                msg = 'pulseAmountIdx of %d is > the number of initial conditions. It must be in the range of initial conditions, 1 to %d.\n';
                error(msgIdent, msg, idx, ode.initialConditionCount);
            end

            idx_previous = idx;
        end
        
    end
    
end % if isfield(ode, 'pulseTimesteps')

if ode.solverTimeLimit <= 0.0
    msgIdent = 'LHS:ODE:SETTINGS:BadSolverTimeLimit';
    msg = 'solverTimeLimit of %f is <= 0.0.\n';
    error(msgIdent, msg, ode.solverTimeLimit);
end

if ode.saveOnError ~= 0 && ode.saveOnError ~= 1
    msgIdent = 'LHS:ODE:SETTINGS:BadSaveOnError';
    msg = 'saveOnError of %d is not 0 or 1.\n';
    error(msgIdent, msg, ode.saveOnError);
end


end




