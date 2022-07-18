function [protocolDescription,channelModel,leakageEC,parameters,solverOptions] = pmBB84Decoy_asymptotic(decoys)
    [protocolDescription,channelModel,leakageEC]=setDescription();
    parameters=setParameters(decoys);
    solverOptions=setOptions();
end

function [protocolDescription,channelModel,leakageEC]=setDescription()
    
    protocolDescription=str2func('pmBB84ReducedLossyDescription');
    channelModel=str2func('pmBB84ReducedDecoyChannel');
    leakageEC=str2func('generalEC');

end

function parameters=setParameters(decoys)

    parameters.names = ["misalignment","loss", "etad","pz","pd","decoys", "f", 'fullstat', 'time', 'ext']; 

	parameters.scan.time = [1];
    parameters.fixed.misalignment = 0;
    parameters.fixed.loss = 0;
    parameters.fixed.etad = 0;
    parameters.fixed.pd = 0;
	parameters.fixed.pz = 0.5;
    parameters.fixed.f = 1;
    parameters.fixed.fullstat = 1;
    parameters.fixed.decoys = decoys;
    parameters.fixed.ext = true;

end

function solverOptions=setOptions()
    solverOptions.globalSetting.cvxSolver = 'mosek';
    solverOptions.globalSetting.cvxPrecision = 'high';
    
    solverOptions.globalSetting.verboseLevel = 1; 
    
    solverOptions.optimizer.name = 'coordinateDescent'; 
    solverOptions.optimizer.linearResolution = 3;
    solverOptions.optimizer.maxIterations = 2; 
    solverOptions.optimizer.linearSearchAlgorithm = 'iterative'; 
    solverOptions.optimizer.iterativeDepth = 2; 
    solverOptions.optimizer.maxSteps = 10; 
    solverOptions.optimizer.optimizerVerboseLevel = 1; 

    solverOptions.solver1.name = 'asymptotic_inequality';
    
    solverOptions.solver1.maxgap = 1e-6; 
    solverOptions.solver1.maxiter = 10;
    solverOptions.solver1.initmethod = 1; 
    
    solverOptions.solver1.linearconstrainttolerance = 1e-10;
    solverOptions.solver1.linesearchprecision = 1e-20;
    solverOptions.solver1.linesearchminstep = 1e-3;
    solverOptions.solver1.maxgap_criteria = true;
    solverOptions.solver1.removeLinearDependence = 'none'; 
    
    solverOptions.solver2.name = 'asymptotic_inequality';
    solverOptions.solver2.epsilon = 0;
    solverOptions.solver2.epsilonprime = 1e-12;
    
    
end

