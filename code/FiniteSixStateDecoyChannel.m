%% FUNCTION NAME: FiniteSixStateDecoyChannel
% Channel model for the 4-6 protocol with decoy, including loss, dark
% count, misalignment, and depolarization.
% Also allows for expectation values to be imported (see the "ext"
% variable).
% Performs finite-size decoy analysis based on Lars' work, which is
% slightly modified from asymptotic decoy analysis

function channelModel = FiniteSixStateDecoyChannel(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["pzB", "pxB", "pzA", "depol", "misalignment", "loss", "etad", "decoys", "pd", ...
        'fullstat', 'time', "eps", "decoyprobs", "ptest", "N"];
    
    %%%%%%%%%%%%%%%%%%%%% interfacing (please do not modify) %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    
    %allocate the outputs of this channel model file
    %can be automatically filled in by calling addExpectations(x) or addExpectations(x,'mask',maskValue)
    expectations = [];
    expMask = [];
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model begin %%%%%%%%%%%%%%%%%%%%%%%%%
    %%fill in calculation process here

    numbases = 2; 
    
    dimA = protocolDescription.dimensions(1);
    dimB = protocolDescription.dimensions(2);
    dimPB = dimB + 1; % dimension of Bob + error

    % axis about which misalignment happens
    axis = [0,1,0];

    % basis choice probabilities for Bob and Alice respectively
    pyB = 1-pzB-pxB;
    bProbB = [pxB, pyB, pzB];
    bProbA = [1-pzA, pzA];
    probList = [ (1-pzA)/2; (1-pzA)/2; pzA/2; pzA/2];

    % get decoy information
%     decoys = [mu1, decoys];
    numDecoys = length(decoys);
%     fprintf("%%%%% OVERRIDING INTENSITY 1 -> %f %%%%%%", mu1)
    
    % Construct entangled state between Alice and Bob in order to determine
    % Alice's expectation values
    phiP = [1;0;0;1]/sqrt(2);
    % apply the loss channel to this ideal state, but with a loss of zero.
    % I do this so that rhoAB has the right dimensions automatically (even
    % if the channel model is changed--as long as the lossChannel function
    % is also changed appropriate :) )
    rhoAB = lossChannel(phiP*phiP', 0);
    observables = protocolDescription.observables;
    % Compute Alice's expectations based on the ideal phiP state to ensure 
    % her state is unchanged
    for iA = 1 : length(observables)
        if protocolDescription.obsMask(iA) == 0
            expectation = trace(observables{iA} * rhoAB);
            addExpectations(expectation)
        end 
    end

%     fprintf('number of expectations before decoy: %d\n', length(expectations))
    
    % simulate the channel (returns 4 x 64 matrix for each decoy)
    %%%%%%%%%%%%%% INSERT EXTERNAL DETECTOR DATA HERE %%%%%%%%%%%%%%%%%%%% 
    dl = DataLoader.instance();
    rawExpectations = dl.getRawExpectations(time);
    
    % convert 64-D pattern to 7-D POVM
    Mapping = createMapping(max(bProbB));
    % mapping is a 64 x 7 matrix

    % bin into squashed model, then do decoy analysis
    bipartiteExpectationsWCP = zeros(4,7,numDecoys);
    for iA = 1:length(decoys)
        bipartiteExpectationsWCP(:,:,iA) = rawExpectations(:,:,iA)*Mapping;
    end
    % bipartiteExpectationsWCP is a 4 x 7 x 3 matrix
        
    LA = size(bipartiteExpectationsWCP, 1);
    LB = size(bipartiteExpectationsWCP, 2);
    % LA should be 4, LB should be 7
    bipartiteExpectationsL = zeros(LA, LB);
    bipartiteExpectationsU = zeros(LA, LB);

    % get solverOptions from caller so we can pass in finite size
    % parameters to decoy analysis
    solverOptions = evalin('caller', 'solverOptions');
    % store all the variables we need for finite decoy analysis in a nice struct
    decoy_parameters.decoyprobs = decoyprobs;
    decoy_parameters.epsPE = eps.PE;
    decoy_parameters.N = N;
    decoy_parameters.options = solverOptions;
    decoy_parameters.ptest = ptest;
    decoy_parameters.bProbA = 0;

    dA = memoize(@decoyAnalysis);

    for iA = 1:LA
        decoy_parameters.bProbA = probList(iA);
        for jB = 1:LB
            decoy_expectations = squeeze(bipartiteExpectationsWCP(iA,jB,:))';
            [bipartiteExpectationsL(iA,jB), bipartiteExpectationsU(iA,jB)] = ...
                dA(decoys, decoy_expectations, decoy_parameters);
        end
    end

    % normalize by signal state probabilities 
    bipartiteExpectationsL = diag(probList) * bipartiteExpectationsL;
    bipartiteExpectationsU = diag(probList) * bipartiteExpectationsU;
    % this will be different if we condition on generation vs test rounds
    
    % QBER and gain statistics
    if(fullstat == 0)
        % constructing coarse grained statistics expectation values
        coarseExpectationsL = zeros(8,1);
        coarseExpectationsU = zeros(8,1);
        itemNum = 0; % iterator to hold position in the above arrays
        for basis = 1:2
            for iA=0:1
                for jB=0:1
                    itemNum = itemNum+1;
                    row = basis+iA+(basis-1);
                    col = basis^2+jB+(basis-1);
                    coarseExpectationsL(itemNum) = bipartiteExpectationsL(row, col);
                    coarseExpectationsU(itemNum) = bipartiteExpectationsU(row, col);
                end
            end
        end
        % normalization (?)
        normL = 1-sum(coarseExpectationsU);
        normU = 1-sum(coarseExpectationsL);
        
        addExpectations(coarseExpectationsL, 'mask', 1);
        addExpectations(normL, 'mask', 1);
        addExpectations(coarseExpectationsU, 'mask', 2);
        addExpectations(normU, 'mask', 2);
    else
        % fine-grained statistics
        fineExpectationsL = zeros(LA*LB, 1);
        fineExpectationsU = zeros(LA*LB, 1);
        itemNum = 0;
        for iRow = 1:LA
            for jCol = 1:LB
                itemNum = itemNum + 1;
                fineExpectationsL(itemNum) = bipartiteExpectationsL(iRow, jCol);
                fineExpectationsU(itemNum) = bipartiteExpectationsU(iRow, jCol);
            end
        end
        addExpectations(fineExpectationsL, 'mask', 1);
        addExpectations(fineExpectationsU, 'mask', 2);
    end
    
    % SIFTING PROBABILITIES FROM DECOY
    signal_simulation = diag(probList) * (rawExpectations(:,:,1)*Mapping);
    % sifting probabilities from decoy: gain = sum of all same-basis
    % events; error is sum of same-basis events that result in different
    % key bits
    gainx  = sum(signal_simulation.*[1,1,0,0,0,0,0; 1,1,0,0,0,0,0; ...
                                     0,0,0,0,0,0,0; 0,0,0,0,0,0,0], 'all');
    errorx = sum(signal_simulation.*[0,1,0,0,0,0,0; 1,0,0,0,0,0,0; ...
                                     0,0,0,0,0,0,0; 0,0,0,0,0,0,0], 'all')/gainx;
    gainz  = sum(signal_simulation.*[0,0,0,0,0,0,0; 0,0,0,0,0,0,0; ...
                                     0,0,0,0,1,1,0; 0,0,0,0,1,1,0], 'all');
    errorz = sum(signal_simulation.*[0,0,0,0,0,0,0; 0,0,0,0,0,0,0; ...
                                     0,0,0,0,0,1,0; 0,0,0,0,1,0,0], 'all')/gainz;

    % signal state proportion (?)
    P1 = decoys(1) * exp(-1*decoys(1));

%     fprintf('number of expectations: %d\n', length(expectations))
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%

    channelModel.expectations = real(expectations);
    channelModel.expMask = expMask;
    channelModel.errorRate = [errorx, errorz];
    channelModel.pSift = [gainx, gainz];
    channelModel.pSignal = P1;

end

%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%

% Functions to define channel operation
% at the moment, depolarizationChannel and rotationChannel are not used.
% They were useful when the protocol was qubit-based.

function rhoPrime = depolarizationChannel(rho, depol) 
    % Nielsen & Chuang p. 379
    krausOps = cell(1,4);
    pauliMatrices = {eye(2), [0,1;1,0], [0,1i;-1i,0], [1,0;0,-1]};
    krausOps{1} = kron(eye(2), sqrt(1-3/4*depol)*eye(2));
    for i = 2 : length(pauliMatrices)
        krausOps{i} = kron(eye(2), sqrt(depol)*pauliMatrices{i}/2);
    end
    rhoPrime = krausFunc(rho, krausOps);
end

function rhoPrime = rotationChannel(rho, theta, axis)
    X = [0,1; 1,0];
    Y = [0, -1i; 1i, 0];
    Z = [1,0; 0,1]; 
    % assume axis is a vector on the bloch sphere
    % make sure it's normalized juuuust in case...
    axis = axis / norm(axis);
    % Rotation by theta about arbitrary axis:
    % cos(theta/2) I - i sin(theta/2) n . sigma
    Ubob = cos(theta)*eye(2) - 1i*sin(theta)*(axis(1)*X + axis(2)*Y + axis(3)*Z);
    % need to include Alice in the kraus operator, even though nothing
    % happens on her system
    U = kron(eye(2), Ubob);
    rhoPrime = krausFunc(rho,{U});
end

function rhoPrime = lossChannel(rho,loss)
    krausOps = cell(1,3);
    % mapping to loss
    krausOps{1} = sqrt(loss)*zket(3,3)*zket(2,1)';
    krausOps{2} = sqrt(loss)*zket(3,3)*zket(2,2)';
    % remainder is not lost
    V = [1,0;0,1;0,0];
    krausOps{3} = sqrt(1-loss)*V;
    %Add alice???s system
    for i = 1:length(krausOps)
        krausOps{i} = kron(eye(2),krausOps{i});
    end
    rhoPrime = krausFunc(rho,krausOps);
end

% Functions to help with decoy analysis

% maps the 64 possible detection events to 7 outcomes:
%               (x+, x-, y+, y-, z+, z-) + null
function mapping = createMapping(px)
    pas = @(prob,N) (1-prob)^N/(2^(N-1)) + prob^N;
    pdis = 1 - 0.5*(pas(px, 3))/(1-pas(px, 3));
    mapping = zeros(64, 7);

    % mapping indices
    x0 = 33; x1 = 17; y0 = 9; y1 = 5; z0 = 3; z1 = 2;
    mappingIndices = [x0, x1, y0, y1, z0, z1];
    sameBasisIndices = [x1, x0, y1, y0, z1, z0];
    % sameBasisIndices(i) is defined such that it is the opposite qubit of 
    % mappingIndices(i) in the same basis ; helpful for double click assignment 

    for iDet = 1 : length(mappingIndices) % which detector
        % single clicks 
        mapping(mappingIndices(iDet), iDet) = 1;
        for oDet = 1 : length(mappingIndices) % the other detector, for cross-clicks
            if oDet ~= iDet
                % calculate the index of the pattern corresponding to both
                % mappingIndices(iDet) and its pair in the same basis;
                % the -1 is because MATLAB is zero indexed :'(
                clickIndex = mappingIndices(iDet) + mappingIndices(oDet) - 1;
                % double clicks -- map to even 0.5/0.5
                if mappingIndices(iDet) == sameBasisIndices(oDet)
                    mapping(clickIndex, iDet) = 0.5;
                else
                    %cross-clicks with two clicks (ex: 101000, 100001, ...)
                    mapping(clickIndex, 7) = pdis;
                    mapping(clickIndex, 1:2) = ones(1,2)*(1-pdis)/2*px;
                    mapping(clickIndex, 3:6) = ones(1,4)*(1-pdis)/2 * (1-px)/2;

                end
            end
        end
    end
    % helper function to count # of 1's in an array
    count1s = @(t) sum(t(:) == 1);
    % fill up the null mapping column and cross-clicks between >2 detectors
    for iPattern = 1 : 64
        % zero detection :)
        if iPattern == 1
            mapping(iPattern, 7) = 1;
        end
        % get the binary representation of this pattern
        pattern = index1to6(iPattern);
        % check if this row corresponds to >2 clicks 
        if count1s(pattern) > 2
            mapping(iPattern, 7) = pdis;
            mapping(iPattern, 1:2) = ones(1,2)*(1-pdis)/2*px;
            mapping(iPattern, 3:6) = ones(1,4)*(1-pdis)/2 * (1-px)/2;
%             mapping(iPattern, 1:6) = ones(1,6)*(1-pdis)/6;
        end
    end
end

% convert a pattern between 1 and 64 to a 1x6 binary array corresponding to
% a detection event
% 2  -> 0,0,0,0,0,1
% 33 -> 1,0,0,0,0,0
function outcome = index1to6(pattern)
    % because i originally wrote this code to work on 
    % numbers from 0 to 63 -_-
    pattern = pattern - 1;
    outcome = zeros(1,6);
    for iBit = 1:6
        e = 6-iBit; % the power of 2 to divide by
        outcome(iBit) = floor(pattern/(2^e));
        pattern = mod(pattern, 2^e);
    end
end

% function to perform decoy state anylsis
function [Y1L, Y1U] = decoyAnalysis(decoys, decoy_expectations, decoy_parameters)
    N = decoy_parameters.N;
    % this is the finite size ball parameter, often referred to as mu
    % however, mu is being used for decoy intensities, so we rename it to lambda here.
    % initial guess for lambda
    lambda0 = sqrt(log(2/decoy_parameters.epsPE)/(2*N));
    % use a rootfinding approximation to make lambda better
    fun = @(x) erfc(-sqrt(2*N)*x) - erfc(sqrt(2*N)*abs(x-1/N)) - 2*(1-eps);
    lambda = fzero(fun,lambda0);
%     fprintf('USING OLD MU FOR NOW')

    % extract the remaining decoy parameters
    t = decoy_parameters.options.solver2.epsilonprime;
    bProbA = decoy_parameters.bProbA;
    decoyprobs = decoy_parameters.decoyprobs;
    ptest = decoy_parameters.ptest;

    cvx_solver mosek

    n_photon = 10; % upper limit of Yn's (# of yields we solve for)
    n_decoy = length(decoys);
    Poisson = @(mu, n) exp(-mu)*mu^n/factorial(n);
    finite_bounds = t + lambda; % when we import data, this will just become lambda I think
    decoy_tolerance = 1e-10;
    Obj = zeros(1, n_photon+1);
    Obj(2) = 1; % Y_1 is the only one we actually want

    % make sure the decoy_expectations are real (or else we have problems!)
    if ~isreal(decoy_expectations)
%         fprintf('**** Decoy expectation values had imaginary components, which are being discarded ****\n');
        decoy_expectations = real(decoy_expectations);
    end

    % solve for upper bound
    try
        cvx_begin quiet
            variable Y(n_photon+1)
            maximize Obj*Y
            for k=1:n_photon + 1
                Y(k) <= 1;
                Y(k) >= 0;
            end
            for i = 1 : n_decoy
                C = zeros(1, n_photon+1); % constraints?
                Ptotal = 0; % running total of the probabilities for the first k photon counts
                for k = 1 : n_photon+1
                    P = Poisson(decoys(i), k-1);
                    C(k) = P;
                    Ptotal = Ptotal + P;
                end
                % p = p(a, mu, test) = probability of sending signal a with
                % intensity choice mu and it's a testing round
                pmu = decoyprobs(i);
                p = bProbA*pmu*ptest;
                % upper equation 10 in decoy paper
                C*Y <= (decoy_expectations(i)*p + finite_bounds)/p + decoy_tolerance;
                C*Y >= (decoy_expectations(i)*p - finite_bounds)/p - (1-Ptotal) - decoy_tolerance;
            end
        cvx_end
        if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
            fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
        end
    catch err
        if(exist('cvx_status', 'var'))
            fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
        else
            fprintf("**** CVX ran into a big issue! ****\n");
            fprintf(err.message);
        end
    end

    Y1U = Y(2);

    % solve for lower bound
    try
        cvx_begin quiet
            variable Y(n_photon + 1)
            minimize Obj * Y
            for k = 1 : n_photon+1
                Y(k) <= 1;
                Y(k) >= 0;
            end
            for i = 1 : n_decoy
                C = zeros(1, n_photon+1);
                Ptotal = 0;
                for k = 1 : n_photon+1
                    P = Poisson(decoys(i), k-1);
                    C(k) = P;
                    Ptotal = Ptotal + P;
                end
                pmu = decoyprobs(i);
                p = bProbA*pmu*ptest;
                C*Y <= (decoy_expectations(i)*p + finite_bounds)/p + decoy_tolerance;
                C*Y >= (decoy_expectations(i)*p - finite_bounds)/p - (1-Ptotal) - decoy_tolerance;
            end
        cvx_end
        if strcmp(cvx_status, 'Infeasible') | strcmp(cvx_status, 'Failed')
            fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
        end
    catch 
        fprintf("**** Warning: decoy state analysis solver exception, status: %s ****\n",cvx_status);
    end

    Y1L = Y(2);
end

