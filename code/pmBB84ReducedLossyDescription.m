%% FUNCTION NAME: SixStateReducedLossyDescription
% Channel description for the 4-6 protocol
%%

function protocolDescription = pmBB84ReducedLossyDescription(names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["pz", 'fullstat'];
    
    %%%%%%%%%%%%%%%%%%%%% interfacing (please do not modify) %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    
    %allocate the outputs of this description file
    %can be automatically filled in by calling addObservables(x) or addObservables(x,'mask',maskValue)
    observables = {};
    obsMask = [];
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description begin %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%fill in calculation process here
    
    dimA = 2;
    dimB = 3;

    numBases = 2; % x and z, sent by Alice

    bProb = [1-pz, pz];

    krausOps = cell(1,2);


    for basis = 1:2%numBases
        % in the qubitProjXYZ function, basis 2 is the Y basis, which we
        % don't want here, so
        if basis == 2
            qubitProjBasis = 3;
        else
            qubitProjBasis = basis; % 1
        end
        %Alice and key register
        tempRA = kron(zket(2,1), qubitProjXYZ(qubitProjBasis, 1, true, 'bra') );
        tempRA = tempRA + kron(zket(2,2), qubitProjXYZ(qubitProjBasis, 2, true, 'bra') );

        % bob and announcements
        
        bob = diag([1,1,0]);
        tempBC = kron(bob, zket(numBases, basis));

        % combine to make kraus operator
        kraus = kron(tempRA, tempBC);

        %add in dependence on basis probability
        kraus = kraus*bProb(basis);
%         disp(size(kraus))
        krausOps{basis} = kraus;
    end

    krausSum = 0;

    for i = 1:length(krausOps)
        krausSum = krausSum + krausOps{i}'*krausOps{i};
    end

    % key projections (i.e. pinching channel)
    keyProj1 = kron(diag([1,0]), eye(dimB*numBases));
    keyProj2 = kron(diag([0,1]), eye(dimB*numBases));
    keyMap={keyProj1, keyProj2};

    % Guarantee Alice's state is unchanged
    basisA = hermitianBasis(dimA);
    for iBasisElm = 1 : length(basisA)
        addObservables(kron(basisA{iBasisElm}, eye(dimB)), 'mask', 0);
    end

    % normalization
    addObservables(eye(dimA*dimB), 'mask', 0);

    h = [1;0];              v = [0;1];
    d = (h+v)/sqrt(2);      a = (h-v)/sqrt(2);
    r = (h+1i*v)/sqrt(2);   l = (h-1i*v)/sqrt(2);
    
    % generate all of Alice's states
    inputStates = {d, a, h, v}; % sending x and z ; change this if different states are sent
    rhoAs = cell(1,4);
    stateChoice = 0; % counter to simplify location in states cell
    for iBasisChoice = 1 : numBases
        for jBitSent = 1 : 2
            stateChoice = stateChoice + 1;
            rhoAs{stateChoice} = bProb(iBasisChoice)*inputStates{stateChoice}*inputStates{stateChoice}';
        end
    end

    % generate all of Bob's states
    % add ;0 for loss dimension
    measureStates = {[d;0], [a;0], [h;0], [v;0]};
    rhoBs = cell(1,5);
    stateChoice = 0;
    for iBasisChoice = 1 : numBases
        for jBitSent = 1 : 2
            stateChoice = stateChoice + 1;
            rhoBs{stateChoice} = bProb(iBasisChoice)*measureStates{stateChoice}*measureStates{stateChoice}';
        end
    end
    % also needs to have a loss dimension POVM element
    rhoBs{end} = [0,0,0;0,0,0;0,0,1];

    if(fullstat==0)
        % add the POVMs that correspond to Alice and Bob having the same basis
        normOp = zeros(6);%normalization observable to be constructed during loop
        for basis = 1:2
            for i=0:1
                for j=0:1
                    row = 2*basis+i-1;
                    col = 2*basis+j-1;
                    bipartiteObs = kron(rhoAs{row}, rhoBs{col});
                    addObservables(bipartiteObs, 'mask', 1);
                    normOp = normOp + bipartiteObs;
                end
            end
        end
    
        % normalization
        addObservables(eye(dimA*dimB)-normOp, 'mask', 1);
    else
        for i = 1 : 4
            for j = 1 : 5
                addObservables(kron(rhoAs{i}, rhoBs{j}), 'mask', 1);
            end
        end
    end

%     fprintf("number of observables: %d\n", length(observables));
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description end %%%%%%%%%%%%%%%%%%%%%%%%%

    protocolDescription.krausOp = krausOps;
    protocolDescription.keyMap = keyMap;
    protocolDescription.observables = observables;
    protocolDescription.obsMask = obsMask;
    protocolDescription.dimensions = [dimA,dimB];

end
