classdef DataLoader < handle
%DataLoader: Class to handle reading expectation value data from an
%external file and dispense it to the program, without having to re-open
%files at each time step. This is a singleton class, so only one instance
%of it ever exists at runtime.
    properties(Access=private)
        file_location
        rawExpectations
        decoys
        parameters
    end
    methods(Access=private)

        % constructor for the DataLoader class. Note that it is private,
        % so it cannot be called directly. Instead, use the
        % DataLoader.instance() method to get an instance
        function obj = DataLoader()
            obj.rawExpectations = [];
            obj.decoys = [];
            obj.parameters = [];
            obj.file_location = '';
        end
    
        % function to actually load data. This will get called any time the
        % file_location field is updated (i.e. it will read the file
        % location whenever is updated.
        function  loadData(obj, file_location)        
            % import the struct of all data
            data = load(file_location);
            % grab the field names to navigate through the data
            fn = fieldnames(data);
        
            decoys_used = zeros(size(fn));
            % grab one experiment from which we grab parameters
            first_experiment = data.(fn{1});
            num_time_steps = length(first_experiment.times);
            mis = first_experiment.misalignment;
            depol = 0; % this isn't in their data structure
            loss = 1-10^(-2.7); % this is also not in the data structure
            pzA = 0.5; % not in the data structure
            pxB = 1 - first_experiment.BS1_transmissivity;
            pzB = first_experiment.BS1_transmissivity * first_experiment.BS2_transmissivity;
            pd = first_experiment.dark_count + first_experiment.background_count;
            obj.parameters = [mis, depol, loss, pzA, pzB, pxB, pd];
        
            % grab the decoy choices so we can pre-allocate the rawExpectations
            % array
            for k = 1 : numel(fn)
                decoys_used(k) = data.(fn{k}).mean_photon_no;
            end
            obj.decoys = unique(decoys_used);
            % discard 0 intensity part for now
            obj.decoys = obj.decoys(obj.decoys > 0);
        
            obj.rawExpectations = zeros(4, 64, length(obj.decoys), num_time_steps);
        
            % for now, grab only the H, V, D, and A (x and z basis) data
            signals_used = ['H','V','D','A'];
            for k = 1 : numel(fn)
                experiment = data.(fn{k});
                % check that this is using the basis data we want
                if(ismember(experiment.signal, signals_used))
                    % discard 0 intensity for now...
                    if(experiment.mean_photon_no > 0)
                        % grab H stuff
                        signal_index = signals_used==experiment.signal;
                        decoy_index = obj.decoys==experiment.mean_photon_no;
                        obj.rawExpectations(signal_index, :, decoy_index, :) = experiment.detections';
                    end
                end
            end
        end
    end
    methods(Static)
        function obj = instance()
            persistent uniqueInstance
            if isempty(uniqueInstance)
                obj = DataLoader();
                uniqueInstance = obj;
            else
                obj = uniqueInstance;
            end
        end
    end
    methods
        function setFileLocation(obj, file_location)
            obj.file_location = file_location;
            obj.loadData(file_location);
        end
        function rExp = getRawExpectations(obj, time)
            if ~isempty(obj.rawExpectations)
                rExp = obj.rawExpectations(:,:,:,time);
            else
                disp('*** ERROR: Attempt to retrieve expectation data when data has not been loaded! ***')
            end
        end
        function decoys = getDecoys(obj)
            decoys = obj.decoys;
        end
        function params = getParameters(obj)
            params = obj.parameters;
        end
    end

end