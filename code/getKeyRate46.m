%% FUNCTION NAME: getKeyRate
% A wrapper around the getKeyRate functionality specifically for the 4-6 protocol. This function 
% assumes the 4-6 procol description, channel model, and error-correction
% leakage, and takes in expectations in the form of raw expectation
% frequencies (before decoy or squashing) and returns upper/lower bound and gap
% This is basically mainIteration.m, but with an assumed protocol
% decsription
%
% Input Data Structure:
% protocol description: [keyMap,krausOperators,observables,obsMask(optional)]
% channel model: [expectations,expMask(optional)probDist/errorRate,pSift]
% EC description: [leakageEC]
% solverOptions: [globalSetting, solver1, solver2]
% (optional) parameter and name list: can be retrieved using findParameter("PARAMETER_NAME",solverOptions)
%
% Output Data Structure:
% lowerBound, upperBound, FWBound, success flag
%%

function result = getKeyRate46(data_location)%[results, parameters]=getKeyRate46(data_location)%
    % Get instance of DataLoader object
    data_obj = DataLoader.instance();
    % Loads the data (could take a bit, esp. for larger data sets)
    data_obj.setFileLocation(data_location);
    % extract decoys and setup parameters from the loaded data
    decoys = data_obj.getDecoys();
    parameters = data_obj.getParameters();
    % set up the preset and assign parameter data
    preset='SixStateDecoy46_asymptotic';
    % parameter extraction starts at 4 because 1-3 are misalignment, loss,
    % and depolarization, which I later realized don't need to be imported.
    % This should be made more sensible in the future.
    etad = parameters(4);   pzA = parameters(5);    pzB = parameters(6);
    pxB = parameters(7);    pd = parameters(8);
    % plug in the extracted decoy and parameter information into the
    % channel description.
    % Note this is quite a bit different from how the software usually
    % assigns parameters--typically, the preset file which contains the
    % input parameters is directly edited, but in this case the parameters
    % are read from the data (or assumed; check the DataLoader class to see
    % and change which are assumed)
    [protocolDescription,channelModel,leakageEC,parameters,solverOptions]=feval(preset, decoys, etad, pzA, pzB, pxB, pd);
    
    helper=helperFunctions; %load helper function file
    [parameters.names,parameters.order]=helper.getOrder(parameters); %generate sorting order or parameters (depending on the given name list)

    %process the scannable (variable) parameters
    parameters_scan = struct2cell(parameters.scan);
    dimensions=helper.getCellDimensions(parameters_scan); %read the dimensions of each parameter to be scanned
    N=prod(dimensions); %total dimension of data to scan

    %process the fixed parameters
    if(~isfield(parameters,'fixed'))
        %no parameters to optimize
        parameters.fixed = struct;
    end
    
    %process the optimizable parameters
    if(~isfield(parameters,'optimize'))
        %no parameters to optimize
        parameters.optimize = struct;
        p_lower=[];
        p_start=[];
        p_upper=[];
        isOptimizing = false;
    else
        optimize_list=cell2mat(struct2cell(parameters.optimize));
        p_lower=optimize_list(:,1)';
        p_start=optimize_list(:,2)';
        p_upper=optimize_list(:,3)';
        isOptimizing = true;
    end

    % %main iteration that test each point in parameters.scan
    % %can optionally uncomment the "parfor" line and comment out "for" to accelerate with multithread/multiple computers
    % %(if parallel computing toolbox is available)
    % %if using parfor, it is advised to set verbose levels to none (here, or in preset file) like below:

%     solverOptions.globalSetting.verboseLevel = 0; %outputs only iteration number
%     solverOptions.optimizer.optimizerVerboseLevel = -1; %outputs nothing when optimizing

%     parfor i=1:N
    for i=1:N
%         if(solverOptions.globalSetting.verboseLevel >= 1)
%             fprintf('main iteration: %d\n',i);
%         end
       
        %%%%%%%%%%%%%% process single-point parameter %%%%%%%%%%%%%%
        indices=helper.expandIndices(i, dimensions); %get indices of parameters.scan
        p_scan = helper.selectCellRow(indices,parameters_scan);
        p_fixed = struct2cell(parameters.fixed)';

        %%%%%% ALSO PROCESS LOSS AS A FUNCTION OF TIME %%%%%%%
        % check if we are scanning over time
        scan_param = fieldnames(parameters.scan);
        scan_param = scan_param{1}; %come on, matlab :(
        if(strcmp(scan_param, 'time'))
            % current time index
            current_time = p_scan;
            % find the index of the loss value in the fixed parameter list
            lossIdx = strcmp(fieldnames(parameters.fixed), 'loss');
            p_fixed(lossIdx) = {[data_obj.getLoss(current_time)]};
        end

        fprintf('time: %d\n',p_scan);
        %%%%%%%%%%%%%% optimize parameter (optional) %%%%%%%%%%%%%%
        
        %optimization of parameters (if parameters.optimize is non-empty)
        if(isOptimizing)
            %wrap the key rate function as a single-input function of "p" (numerical array of optimizable parameters)
            rateFunction = @(p) getKeyRate_wrapper(parameters.names,helper.reorder([p_scan,p_fixed,p],parameters.order),protocolDescription,channelModel,leakageEC,solverOptions); %convert getKeyRate as a single-value function f(p) for optimization
            if(solverOptions.globalSetting.verboseLevel >= 1)
                fprintf('begin optimization\n');
            end
            if(strcmp(solverOptions.optimizer.name,'bruteForce'))
                p_optimal = helper.bruteForceSearch(rateFunction,p_lower,p_upper,solverOptions.optimizer);
            elseif(strcmp(solverOptions.optimizer.name,'coordinateDescent'))
                p_optimal = helper.coordinateDescent(rateFunction,p_start,p_lower,p_upper,solverOptions.optimizer);
            elseif(strcmp(solverOptions.optimizer.name,'gradientDescent'))
                p_optimal = helper.gradientDescent(rateFunction,p_start,p_lower,p_upper,solverOptions.optimizer);
            elseif(strcmp(solverOptions.optimizer.name,'localSearch_Adam'))
                p_optimal = helper.localSearch_Adam(rateFunction,p_start,p_lower,p_upper,solverOptions.optimizer);
            end
            if(solverOptions.globalSetting.verboseLevel >= 1)
                fprintf('finished optimization\n');
            end
        else
            p_optimal = [];
        end

        %%%%%%%%%%%%%% evaluate descriptions %%%%%%%%%%%%%%
        
        %generation of single-row parameter list (a cell array)
        p_scan = num2cell(p_scan); %convert p_scan to cell array
        if(~isempty(p_optimal))
            p_optimal = num2cell(p_optimal); %convert p_optimal to cell array
        end
        p_full=[p_scan,p_fixed,p_optimal]; %concatenate with p_fixed
        p_full = helper.reorder(p_full,parameters.order);
        
        %%%%%%%%%%%%%% optimize loss, pxB, and pzB based on simulation data
        %%%%% after some tests, it looks like this doesn't help key rate
        %%%%% much and only adds probabilistic instability
        % only if we are scanning over time!!
%         scan_params = fieldnames(parameters.scan);
%         if(scan_params{1} == 'time')
%             lossIdx = parameters.names=='loss';
%             pxBIdx = parameters.names=='pxB';
%             pzBIdx = parameters.names=='pzB';
%             
%             % p_scan is the current time
%             [optLoss, optpzB, optpxB] = optimizeParameters(decoys, mis, depol, loss, etad, pzA, pzB, pxB, pd, p_scan{1});
%             p_full{lossIdx} = optLoss;
%             p_full{pxBIdx} = optpxB;
%             p_full{pzBIdx} = optpzB;
%         end
        %%%%%%%%%%%%%% end parameter optimization

        %evaluate the protocol description, channel model, and leakage
        thisProtocolDescription=protocolDescription(parameters.names,p_full);
        thisChannelModel=channelModel(thisProtocolDescription,parameters.names,p_full);
        thisLeakageEC=leakageEC(thisChannelModel,parameters.names,p_full);
        

        %%%%%%%%%%%%%% perform calculation %%%%%%%%%%%%%%
        
        %calculate key rate by calling the solver module
        %note that the full parameter list is also passed into the function for the solver to optionally directly access (e.g. security parameter eps, data size N, etc.).
%         [results(i).lowerBound,results(i).upperBound,results(i).FWBound,results(i).debugInfo] = getKeyRate(thisProtocolDescription,thisChannelModel,thisLeakageEC,solverOptions,p_full,parameters.names);
         [results{i}.lowerBound,results{i}.upperBound,results{i}.FWBound,results{i}.debugInfo] = getKeyRate(thisProtocolDescription,thisChannelModel,thisLeakageEC,solverOptions,p_full,parameters.names);

        %also save the current parameter set for reference (1-D cell array)
%         results(i).debugInfo.current_parameters = p_full;
%         results(i).debugInfo.names = parameters.names;
        results{i}.debugInfo.current_parameters = p_full;
        results{i}.debugInfo.names = parameters.names;
    end
    
    result.results = results;
    result.parameters = parameters;
end


%%% function to optimize the parameters based on imported data; in
%%% particular, this is meant to make Alice's part of the system line up
%%% with the parameters of the experiment found in the raw detector data.
% simulates the coherent source channel for each decoy for the
% optimization procedure.
function [loss, pz, px] = optimizeParameters(decoys, misalignment, depol, loss, etad, pzA, pzB, pxB, pd, time)
    dl = DataLoader.instance();
    imported_expectations = dl.getRawExpectations(time);
    x0 = [loss, pzB, pxB];
    a = {decoys, etad, misalignment, [0,1,0], depol, pd};
    % define the function we want to minimize; a is a cell array with
    % parameters that aren't changed (decoys, detector efficiency, dark
    % count rate, etc.) and x is a regular array with parameters that CAN
    % be changed (loss, depolarization parameter, basis choice
    % probabilities for Bob)
    f = @(x,a) minimizeDistance(simChannel(a{1}, x(1), a{2}, a{3}, a{4}, a{5}, x(2), x(3), a{6}), imported_expectations);
    % now define the actual optimization function, which is a function of x
    % only
    optF = @(x) f(x,a);
    % perform the search
    x = fminsearch(optF, x0);
    % update parameters if they make sense
    if(x(1) < 1 && x(1) >= 0)
        loss = x(1);
    end
    if(x(2) < 1 && x(2) > 0 && x(3) < 1 && x(3) > 0)
        pz = x(2);
        px = x(3);
    else
        pz = pzB;
        px = pxB;
    end
end

%%%%%%% helper functions for parameter optimization %%%%%%%%%%
function sim = simChannel(decoys, loss, etad, mis, axis, depol, pzB, pxB, pd)
    sim = zeros(4,64,length(decoys));
    for i = 1 : length(decoys)
        sim(:,:,i) = coherentSourceChannel(decoys(i), loss, etad, mis, axis, depol, pzB, pxB, pd);
    end
end

% objective function for the optimization procedure to fit to experimental
% parameters; currently finds the absolute value of the maximum difference,
% but this may not be the best distance measure.
function objective = minimizeDistance(data1, data2)
%     obj = 0;
%     for i = 1 : size(data1,3)
%         obj = obj + norm(data1(:,:,i)-data2(:,:,i));
%     end
%     objective = obj;
    objective = max(max(max(abs(data1-data2))));
end
