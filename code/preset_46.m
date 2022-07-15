
function [protocolDescription,channelModel,leakageEC,parameters,solverOptions] = SixStateDecoy46_asymptotic(decoys)
    [protocolDescription,channelModel,leakageEC]=setDescription();
    parameters=setParameters(decoys);
    solverOptions=setOptions();
end
function [protocolDescription,channelModel,leakageEC]=setDescription()
    
    protocolDescription=str2func('SixStateReducedLossyDescription');
    channelModel=str2func('SixStateDecoyChannel'); 
    leakageEC=str2func('generalEC');

end
function parameters=setParameters(decoys)

    parameters.names = ["etad", "pzB","pzA", "pxB","pd","decoys", "f", 'fullstat', 'time']; 
    parameters.scan.time = [342 343 344 345 346 347 348 349 350 351 352];
    parameters.fixed.pzA = 0.5; 
    parameters.fixed.pzB = 0.167;
    parameters.fixed.pxB = 0.666;
    parameters.fixed.pd = 1e-6;
    parameters.fixed.f = 1.16;
    parameters.fixed.fullstat = 1;
    parameters.fixed.etad = 1;
    parameters.fixed.decoys = decoys;
end

function solverOptions=setOptions()
    solverOptions.globalSetting.cvxSolver = 'mosek';
    solverOptions.globalSetting.cvxPrecision = 'high';
    
    solverOptions.globalSetting.verboseLevel = 2; 
  
    solverOptions.optimizer.name = 'coordinateDescent'; 
    solverOptions.optimizer.linearResolution = 3; 
    solverOptions.optimizer.maxIterations = 2; 
    solverOptions.optimizer.linearSearchAlgorithm = 'iterative'; 
    solverOptions.optimizer.iterativeDepth = 2; 
    solverOptions.optimizer.maxSteps = 10; 
    solverOptions.optimizer.optimizerVerboseLevel = 0; 

    solverOptions.solver1.name = 'asymptotic_inequality';
    solverOptions.solver1.maxgap = 1e-6; 
    solverOptions.solver1.maxiter = 20;%10;
    solverOptions.solver1.initmethod = 1; %minimizes norm(rho0-rho) or -lambda_min(rho), use method 1 for finite size, 2 for asymptotic v1
    solverOptions.solver1.linearconstrainttolerance = 1e-8;
    solverOptions.solver1.linesearchprecision = 1e-20;
    solverOptions.solver1.linesearchminstep = 1e-3;
    solverOptions.solver1.maxgap_criteria = false; %true for testing gap, false for testing f1-f0
    solverOptions.solver1.removeLinearDependence = 'none'; %'qr' for QR decomposition, 'rref' for row echelon form, empty '' for not removing linear dependence
    solverOptions.solver2.name = 'asymptotic_inequality';
    solverOptions.solver2.epsilon = 0;
    solverOptions.solver2.epsilonprime = 1e-12;
end
