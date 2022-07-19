function [protocolDescription,channelModel,leakageEC,parameters,solverOptions] = pm46Decoy_finite(decoys, N)
    [protocolDescription,channelModel,leakageEC]=setDescription();
    parameters=setParameters(decoys, N);
    solverOptions=setOptions();
end

function [protocolDescription,channelModel,leakageEC]=setDescription()
    
    protocolDescription=str2func('SixStateReducedLossyDescription');
    channelModel=str2func('FiniteSixStateDecoyChannel');
    leakageEC=str2func('generalEC');

end

function parameters=setParameters(decoys, N)

    parameters.names = ["pzB","pzA", "pxB","decoys", "f", ...
        'fullstat', 'time', "N","ptest","eps","alphabet","postselection","physDimAB", "decoyprobs"]; %BB84 Decoy

    parameters.scan.time = 340:350;
    
    parameters.fixed.pzA = 0.5; %basis choice probability (for Z basis)
    % we found that choosing the following values for pzB and pxB works
    % much better for getting key rate than the original parameters of 
    % px = 0.5, pz = 0.25
    parameters.fixed.pzB = 1/6;
    parameters.fixed.pxB = 2/3;
    parameters.fixed.f = 1.16;
    parameters.fixed.fullstat = 1;
    parameters.fixed.decoys = decoys; %first is signal, other two are decoys
    parameters.fixed.decoyprobs = ones(size(decoys))/length(decoys);

    % parameters for finite-size analysis
    parameters.fixed.N = N; %(finite size) sent data
    parameters.fixed.ptest = parameters.fixed.N^(-1/8); %(finite size) proportion of data used in testing
    parameters.fixed.alphabet = 2; %(finite size) the encoding alphabet size - for qubits the size is 2

    eps.PE = (1/4)*1e-8; %failure probability for parameter estimation 
    eps.bar = (1/4)*1e-8; %uncertainty from "smoothing" of min-entropy among states similar to feasible state
    eps.EC = (1/4)*1e-8; %failure probability for error-correction
    eps.PA = (1/4)*1e-8; %failure probability for privacy amplification
    parameters.fixed.eps = eps;

    %settings for optional post-selection technique in finite-size analysis
    %note that usually eps has to be very small (e.g. 1e-200) in order for post-selection to return key rate
    parameters.fixed.postselection = 0; %(finite size) whether to use or not use post-selection, corresponding to coherent/collective attack
    parameters.fixed.physDimAB = 2*3; %(finite size) physical dimensions of Alice's and Bob's outgoing/incoming signals, used for post-selection

end

function solverOptions=setOptions()
    %%%%%%%%%%%%%%%%% global setting %%%%%%%%%%%%%%%%%
    solverOptions.globalSetting.cvxSolver = 'mosek'; 
    solverOptions.globalSetting.cvxPrecision = 'medium';
    
    %output level:
    %0.output nothing (except errors)
    %1.output at each main iteration
    %2.output at each solver 1 FW iteration
    %3.show all SDP solver output
    solverOptions.globalSetting.verboseLevel = 1; 
    
    % key rate calculation selection
    % empty to default to whatever the solver is
    % 'asymptotic' for no finite size corrections
    % 'finite' for finite size corrections
    solverOptions.globalSetting.keyRateFormula = 'finite';

    %%%%%%%%%%%%%%%%% parameter optimizer setting %%%%%%%%%%%%%%%%%
    solverOptions.optimizer.name = 'coordinateDescent'; %choose between 'coordinateDescent' and 'bruteForce'
    solverOptions.optimizer.linearResolution = 20; %resolution in each dimension (for brute force search and coordinate descent)
    solverOptions.optimizer.maxIterations = 3; %max number of iterations (only for coordinate descent)
    solverOptions.optimizer.linearSearchAlgorithm = 'iterative'; %choose between built-in 'fminbnd' and custom 'iterative' algorithms for linear search (only for coordinate descent)
    solverOptions.optimizer.iterativeDepth = 2; %choose depth of iteration levels; function count = depth * linearResolution (only for coordinate descent and if using 'iterative')
    solverOptions.optimizer.maxSteps = 20; %max number of steps (only for gradient descent and ADAM)
    solverOptions.optimizer.optimizerVerboseLevel = 1; %0:only output optimized result; 1:output a progress bar; 2:output at each function call

    %%%%%%%%%%%%%%%%% step1Solver setting %%%%%%%%%%%%%%%%%
    solverOptions.solver1.name = 'asymptotic_inequality';
    
    %options mainly affecting performance
    solverOptions.solver1.maxgap = 2.5e-3;%1e-6; %1e-6 for asymptotic, 2.5e-3 for finite;
    solverOptions.solver1.maxiter = 10;
    solverOptions.solver1.initmethod = 1; %minimizes norm(rho0-rho) or -lambda_min(rho), use method 1 for finite size, 2 for asymptotic v1
    
    %default options
    solverOptions.solver1.linearconstrainttolerance = 1e-10;
    solverOptions.solver1.linesearchprecision = 1e-20;
    solverOptions.solver1.linesearchminstep = 1e-3;
    solverOptions.solver1.maxgap_criteria = false; %true for testing gap, false for testing f1-f0
    solverOptions.solver1.removeLinearDependence = 'rref'; %'qr' for QR decomposition, 'rref' for row echelon form, empty '' for not removing linear dependence

    %%%%%%%%%%%%%%%%% step2Solver setting %%%%%%%%%%%%%%%%%
    solverOptions.solver2.name = 'asymptotic_inequality';
    solverOptions.solver2.epsilon = 0;
    solverOptions.solver2.epsilonprime = 1e-10;
    solverOptions.solver2.linearconstrainttolerance = solverOptions.solver1.linearconstrainttolerance;
    
end

