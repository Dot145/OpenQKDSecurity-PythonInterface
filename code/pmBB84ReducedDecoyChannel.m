%% FUNCTION NAME: SixStateDecoyChannel
% Channel model for the 4-6 protocol with decoy, including loss, dark
% count, misalignment, and depolarization.
%%

function channelModel = pmBB84ReducedDecoyChannel(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["pz", "misalignment", "loss", "etad", "decoys", "pd", 'fullstat', 'time', 'ext'];
    
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

    % basis choice probabilities for Bob and Alice respectively
    px = 1-pz;
    bProb = [1-pz, px];
    probList = [ (1-pz)/2; (1-pz)/2; pz/2; pz/2];

    % get decoy information
    numDecoys = length(decoys);
    
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
    for i = 1 : length(observables)
        if protocolDescription.obsMask(i) == 0
            expectation = trace(observables{i} * rhoAB);
            addExpectations(expectation)
        end 
    end

%     fprintf('number of expectations before decoy: %d\n', length(expectations))
    
    % simulate the channel (returns 4 x 64 matrix for each decoy)
    %%%%%%%%%%%%%% INSERT EXTERNAL DETECTOR DATA HERE %%%%%%%%%%%%%%%%%%%% 
    dl = DataLoader.instance();
    rawExpectations = dl.getRawExpectations(time);
    
    % convert 16-D pattern to 5-D POVM
    Mapping = createMapping();
    % mapping is a 16 x 5 matrix

    % bin into squashed model, then do decoy analysis
    bipartiteExpectationsWCP = zeros(4,5,numDecoys);
    for i = 1:length(decoys)
        bipartiteExpectationsWCP(:,:,i) = rawExpectations(:,:,i)*Mapping;
    end
    % bipartiteExpectationsWCP is a 4 x 5 x 3 matrix
        
    L1 = size(bipartiteExpectationsWCP, 1);
    L2 = size(bipartiteExpectationsWCP, 2);
    % i think L1 should be 4, L2 should be 5
    bipartiteExpectationsL = zeros(L1, L2);
    bipartiteExpectationsU = zeros(L1, L2);

    dA = memoize(@decoyAnalysis);

    for i = 1:L1
        for j = 1:L2
            decoy_expectations = squeeze(bipartiteExpectationsWCP(i,j,:))';
            [bipartiteExpectationsL(i,j), bipartiteExpectationsU(i,j)] = ...
                dA(decoys, decoy_expectations);
        end
    end

    % normalize by signal state probabilities (?)
    bipartiteExpectationsL = diag(probList) * bipartiteExpectationsL;
    bipartiteExpectationsU = diag(probList) * bipartiteExpectationsU;
    
    % QBER and gain statistics (?)
    if(fullstat == 0)
        % constructing coarse grained statistics expectation values
        coarseExpectationsL = zeros(8,1);
        coarseExpectationsU = zeros(8,1);
        itemNum = 0; % iterator to hold position in the above arrays
        for basis = 1:2
            for i=0:1
                for j=0:1
                    itemNum = itemNum+1;
                    row = basis+i+(basis-1);
                    col = basis^2+j+(basis-1);
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
        fineExpectationsL = zeros(L1*L2, 1);
        fineExpectationsU = zeros(L1*L2, 1);
        itemNum = 0;
        for iRow = 1:L1
            for jCol = 1:L2
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
    % 
    gainx  = sum(signal_simulation.*[1,1,0,0,0; 1,1,0,0,0; ...
                                     0,0,0,0,0; 0,0,0,0,0], 'all');
    errorx = sum(signal_simulation.*[0,1,0,0,0; 1,0,0,0,0; ...
                                     0,0,0,0,0; 0,0,0,0,0], 'all')/gainx;
    gainz  = sum(signal_simulation.*[0,0,0,0,0; 0,0,0,0,0; ...
                                     0,0,1,1,0; 0,0,1,1,0], 'all');
    errorz = sum(signal_simulation.*[0,0,0,0,0; 0,0,0,0,0; ...
                                     0,0,0,1,0; 0,0,1,0,0], 'all')/gainz;

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
    %Add alice’s system
    for i = 1:length(krausOps)
        krausOps{i} = kron(eye(2),krausOps{i});
    end
    rhoPrime = krausFunc(rho,krausOps);
end

% Functions to help with decoy analysis

% maps the 16 possible detection events to 5 outcomes:
%               (x+, x-, z+, z-) + null
function mapping = createMapping()
    %Patterns:
    %column: [H,V,+,-,None]
    %row: Alice POVM
    MappingH=zeros(16,1);
    MappingH(9)=1; %1000
    MappingH(13)=0.5; %1100
    MappingV=zeros(16,1);
    MappingV(5)=1; %0100
    MappingV(13)=0.5; %1100
    
    MappingD=zeros(16,1);
    MappingD(3)=1; %0010
    MappingD(4)=0.5; %0011
    
    MappingA=zeros(16,1);
    MappingA(2)=1; %0001
    MappingA(4)=0.5; %0011
    
    MappingN=zeros(16,1);
    MappingN(1)=1; %0000
    MappingN(6)=1; %0101
    MappingN(7)=1; %0110
    MappingN(8)=1; %0111
    MappingN(10)=1; %1001
    MappingN(11)=1; %1010
    MappingN(12)=1; %1011
    MappingN(14)=1; %1101
    MappingN(15)=1; %1110
    MappingN(16)=1; %1111
    
    mapping=[MappingH,MappingV,MappingD,MappingA,MappingN];
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
% keeping this bit just in case everything breaks, but it has weird
% ordering (000000 -> 100000 -> 010000 ...)
% function outcome = index1to6(pattern)
%     % because i originally wrote this code to work on 
%     % numbers from 0 to 63 -_-
%     pattern = pattern - 1;
%     outcome = zeros(1,6);
%     for iBit = 1:6
%         e = 6-iBit; % the power of 2 to divide by
%         outcome(e+1) = floor(pattern/(2^e));
%         pattern = mod(pattern, 2^e);
%     end
% end

% function to perform decoy state anylsis
function [Y1L, Y1U] = decoyAnalysis(decoys, decoy_expectations)
    cvx_solver mosek

    n_photon = 10; % upper limit of Yn's (# of yields we solve for)
    n_decoy = length(decoys);
    Poisson = @(mu, n) exp(-mu)*mu^n/factorial(n);
    decoy_tolerance = 1e-12;
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
            
                % upper equation 10 in decoy paper
                C*Y <= decoy_expectations(i) + decoy_tolerance;
                C*Y >= decoy_expectations(i) - decoy_tolerance - (1-Ptotal);
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
                C*Y <= decoy_expectations(i) + decoy_tolerance;
                C*Y >= decoy_expectations(i) - decoy_tolerance - (1-Ptotal);
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

%helper function that simulates the channel for WCP sources for BB84
function expectations=coherentSourceChannelBB84(mu,eta,etad,pd,ed,px)
    expectations=zeros(4,16);
    pz=1-px;
    t=eta*etad; %total transmittance
    Poisson=@(mu,n) exp(-mu)*mu^n/factorial(n);

    theta=ed;
    PL=sin(pi/4-theta)^2;
    PU=cos(pi/4-theta)^2;
    % see page 9 of the paper
    mapping_passive=[[pz*(1-ed),pz*ed,px*PU,px*PL];[pz*ed,pz*(1-ed),px*PL,px*PU];[pz*PL,pz*PU,px*(1-ed),px*ed];[pz*PU,pz*PL,px*ed,px*(1-ed)]];
    
    for input=1:4
        %iterating over each input state
        
        for output=1:16
            %iterating over each pattern
            a=index1to4(output-1); %detector event, 4 elements corresponding to [H,V,D,A]
            
            Ppattern=1;
            %passive basis choice
            for k=1:4
                %iterating over each detector
                Pclick=1-Poisson(mu*t*mapping_passive(input,k),0); 
                % effect of dark count
                Pclick=1-(1-Pclick)*(1-pd); 
                if(a(k)==1)
                    Ppattern=Ppattern*Pclick;
                elseif(a(k)==0)
                    Ppattern=Ppattern*(1-Pclick);
                end
            end
            expectations(input,output)=Ppattern;
        end
    end
end

%takes in an index of 1 and convert to array of 4
function a=index1to4(index)
    a(1) = floor(index/8);
    index = mod(index,8);
    a(2) = floor(index/4);
    index = mod(index,4);
    a(3) = floor(index/2);
    index = mod(index,2);
    a(4) = index;
end
