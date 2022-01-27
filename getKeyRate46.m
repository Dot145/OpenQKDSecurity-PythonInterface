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
    mis = parameters(1);    depol = parameters(2);  loss = parameters(3);
    pzA = parameters(4);    pzB = parameters(5);    pxB = parameters(6);
    pd = parameters(7);
    % plug in the extracted decoy and parameter information into the
    % channel description.
    % Note this is quite a bit different from how the software usually
    % assigns parameters--typically, the preset file which contains the
    % input parameters is directly edited, but in this case the parameters
    % are read from the data (or assumed; check the DataLoader class to see
    % and change which are assumed)
    [protocolDescription,channelModel,leakageEC,parameters,solverOptions]=feval(preset, decoys, mis, depol, loss, pzA, pzB, pxB, pd);
    
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
    disp(size(result))
end


% old data loading functions; this functionality was replaced by DataLoader
% function [rawExpectations, decoys, parameters] = loadData()
%     num_time_steps = 0;
% 
%     % import the struct of all data
%     data = load('data/Waterloo_fullData/alldata.mat');
%     % grab the field names to navigate through the data
%     fn = fieldnames(data);
% 
%     decoys_used = zeros(size(fn));
%     % grab one experiment from which we grab parameters
%     first_experiment = data.(fn{1});
%     num_time_steps = length(first_experiment.times);
%     mis = first_experiment.misalignment;
%     depol = 0; % this isn't in their data structure
%     loss = 1-10^(-2.7); % this is also not in the data structure
%     pzA = 0.5; % not in the data structure
%     pxB = 1 - first_experiment.BS1_transmissivity;
%     pzB = first_experiment.BS1_transmissivity * first_experiment.BS2_transmissivity;
%     pd = first_experiment.dark_count + first_experiment.background_count;
%     parameters = [mis, depol, loss, pzA, pzB, pxB, pd];
% 
%     % grab the decoy choices so we can pre-allocate the rawExpectations
%     % array
%     for k = 1 : numel(fn)
%         decoys_used(k) = data.(fn{k}).mean_photon_no;
%     end
%     decoys = unique(decoys_used);
%     % discard 0 intensity part for now
%     decoys = decoys(decoys > 0);
% 
%     rawExpectations = zeros(4, 64, length(decoys), num_time_steps);
% 
%     % for now, grab only the H, V, D, and A (x and z basis) data
%     signals_used = ['H','V','D','A'];
%     for k = 1 : numel(fn)
%         experiment = data.(fn{k});
%         % check that this is using the basis data we want
%         if(ismember(experiment.signal, signals_used))
%             % discard 0 intensity for now...
%             if(experiment.mean_photon_no > 0)
%                 % grab H stuff
%                 signal_index = signals_used==experiment.signal;
%                 decoy_index = decoys==experiment.mean_photon_no;
%                 rawExpectations(signal_index, :, decoy_index, :) = experiment.detections';
%             end
%         end
%     end
% end
% 
% function rawExpectations = loadDataOld(num_decoys, time)
%     rawExpectations = zeros(4, 64, num_decoys);
%     money_row = time + 347; % the column of the detections matrix that corresponds to time = 0 is 347
%     % decoy 1: mu = 0.5
%     mu5dat_x0 = load('data/Waterloo_fullData/Zenith_27dB_rx1ry0rz0_MeanPhotonNo0.5_ClickProb.mat');
%     rawExpectations(1,:,1) = mu5dat_x0.detections(money_row,:);
%     mu5dat_x1 = load('data/Waterloo_fullData/Zenith_27dB_rx-1ry0rz0_MeanPhotonNo0.5_ClickProb.mat');
%     rawExpectations(2,:,1) = mu5dat_x1.detections(money_row,:);
%     mu5dat_z0 = load('data/Waterloo_fullData/Zenith_27dB_rx0ry0rz1_MeanPhotonNo0.5_ClickProb.mat');
%     rawExpectations(3,:,1) = mu5dat_z0.detections(money_row,:);
%     mu5dat_z1 = load('data/Waterloo_fullData/Zenith_27dB_rx0ry0rz-1_MeanPhotonNo0.5_ClickProb.mat');
%     rawExpectations(4,:,1) = mu5dat_z1.detections(money_row,:);
%     % decoy 2: mu = 0.7
%     mu7dat_x0 = load('data/Waterloo_fullData/Zenith_27dB_rx1ry0rz0_MeanPhotonNo0.7_ClickProb.mat');
%     rawExpectations(1,:,2) = mu7dat_x0.detections(money_row,:);
%     mu7dat_x1 = load('data/Waterloo_fullData/Zenith_27dB_rx-1ry0rz0_MeanPhotonNo0.7_ClickProb.mat');
%     rawExpectations(2,:,2) = mu7dat_x1.detections(money_row,:);
%     mu7dat_z0 = load('data/Waterloo_fullData/Zenith_27dB_rx0ry0rz1_MeanPhotonNo0.7_ClickProb.mat');
%     rawExpectations(3,:,2) = mu7dat_z0.detections(money_row,:);
%     mu7dat_z1 = load('data/Waterloo_fullData/Zenith_27dB_rx0ry0rz-1_MeanPhotonNo0.7_ClickProb.mat');
%     rawExpectations(4,:,2) = mu7dat_z1.detections(money_row,:);
%     % decoy 3: mu = 0.9
%     mu9dat_x0 = load('data/Waterloo_fullData/Zenith_27dB_rx1ry0rz0_MeanPhotonNo0.9_ClickProb.mat');
%     rawExpectations(1,:,3) = mu9dat_x0.detections(money_row,:);
%     mu9dat_x1 = load('data/Waterloo_fullData/Zenith_27dB_rx-1ry0rz0_MeanPhotonNo0.9_ClickProb.mat');
%     rawExpectations(2,:,3) = mu9dat_x1.detections(money_row,:);
%     mu9dat_z0 = load('data/Waterloo_fullData/Zenith_27dB_rx0ry0rz1_MeanPhotonNo0.9_ClickProb.mat');
%     rawExpectations(3,:,3) = mu9dat_z0.detections(money_row,:);
%     mu9dat_z1 = load('data/Waterloo_fullData/Zenith_27dB_rx0ry0rz-1_MeanPhotonNo0.9_ClickProb.mat');
%     rawExpectations(4,:,3) = mu9dat_z1.detections(money_row,:);
%     if(num_decoys > 3)
%         % decoy 4: mu = 0.
%         mu0dat_x0 = load('data/Waterloo_fullData/Zenith_27dB_rx1ry0rz0_MeanPhotonNo0._ClickProb.mat');
%         rawExpectations(1,:,4) = mu0dat_x0.detections(money_row,:);
%         mu0dat_x1 = load('data/Waterloo_fullData/Zenith_27dB_rx-1ry0rz0_MeanPhotonNo0._ClickProb.mat');
%         rawExpectations(2,:,4) = mu0dat_x1.detections(money_row,:);
%         mu0dat_z0 = load('data/Waterloo_fullData/Zenith_27dB_rx0ry0rz1_MeanPhotonNo0_ClickProb.mat');
%         rawExpectations(3,:,4) = mu0dat_z0.detections(money_row,:);
%         mu0dat_z1 = load('data/Waterloo_fullData/Zenith_27dB_rx0ry0rz-1_MeanPhotonNo0_ClickProb.mat');
%         rawExpectations(4,:,4) = mu0dat_z1.detections(money_row,:);
%     end
% end
