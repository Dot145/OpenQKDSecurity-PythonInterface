
function [protocolDescription,channelModel,leakageEC,parameters,solverOptions] = SixStateDecoy46_asymptotic(decoys, mis, depol, loss, etad, pzA, pzB, pxB, pd)
    [protocolDescription,channelModel,leakageEC]=setDescription();
    parameters=setParameters(decoys, mis, depol, loss, etad, pzA, pzB, pxB, pd);
    solverOptions=setOptions();
end
function [protocolDescription,channelModel,leakageEC]=setDescription()
    
    protocolDescription=str2func('SixStateReducedLossyDescription');
    channelModel=str2func('SixStateDecoyChannel'); 
    leakageEC=str2func('generalEC');

end
function parameters=setParameters(decoys, mis, depol, loss, etad, pzA, pzB, pxB, pd)

    parameters.names = ["misalignment","loss", "etad", "depol","pzB","pzA", "pxB","pd","decoys", "f", 'fullstat', 'time', 'ext']; %BB84 Decoy
    parameters.scan.time = [342 344 346 348 350 352];
    parameters.fixed.misalignment = mis;
    parameters.fixed.depol = depol;
    parameters.fixed.pzA = pzA; 
    % we found that choosing the following values for pzB and pxB works
    % much better for getting key rate than the original parameters of 
    % px = 0.5, pz = 0.25
    parameters.fixed.pzB = 0.167;
    parameters.fixed.pxB = 0.666;
    parameters.fixed.pd = pd;
    parameters.fixed.f = 1;
    parameters.fixed.fullstat = 1;
    parameters.fixed.loss = loss
    parameters.fixed.etad = etad;
    parameters.fixed.decoys = decoys;
    parameters.fixed.ext = true;    
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
    solverOptions.optimizer.optimizerVerboseLevel = 1; 

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
