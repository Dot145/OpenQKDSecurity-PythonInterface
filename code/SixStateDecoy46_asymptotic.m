%Preset input for prepare-and-measure BB84 with asymptotic solver
%There is squashing model for Bob's POVM, and channel model contains loss, misalignment, and dark count.
%We assume Alice uses WCP source (where only single-photon portion is used to generate key)
%A decoy-state analysis is performed in the channel model file after generating expectations, to estimate single-photon contributions
%In this preset we are using the asymptotic (decoy version) solver on the 4-6 protocol.

function [protocolDescription,channelModel,leakageEC,parameters,solverOptions] = SixStateDecoy46_asymptotic(decoys, mis, depol, loss, pzA, pzB, pxB, pd)
    [protocolDescription,channelModel,leakageEC]=setDescription();
    parameters=setParameters(decoys, mis, depol, loss, pzA, pzB, pxB, pd);
    solverOptions=setOptions();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% user input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set description files
%returns function handle protocolDescription(parameters) and expectations(parameters); description is a struct with multiple fields
%string name must match file/function name; case-sensitive
function [protocolDescription,channelModel,leakageEC]=setDescription()
    
    protocolDescription=str2func('SixStateReducedLossyDescription');
    channelModel=str2func('SixStateDecoyChannel'); %contains asymptotic decoy state analysis in file
    leakageEC=str2func('generalEC');

end

%set the input parameters
%three types of parameters are considered: scan, fixed, optimize
function parameters=setParameters(decoys, mis, depol, loss, pzA, pzB, pxB, pd)

    parameters.names = ["misalignment","loss","depol","pzB","pzA", "pxB","pd","decoys", "f", 'fullstat', 'time', 'ext']; %BB84 Decoy

    %%%%%%%%%%%%%%%% 1.parameter to scan over %%%%%%%%%%%%%%%%
    %must name at least one parameter to scan (can be a single-point array if only interested in a fixed value)
    
    parameters.scan.time =  (-2:2) + 347;%t=0 is row 347

    %%%%%%%%%%%%%%%% 2.fixed parameters %%%%%%%%%%%%%%%%
    %optional; the constant values can be either numerical or array/matrices
    
    parameters.fixed.misalignment = mis;%0.01; 
    parameters.fixed.depol = depol;%0.01;
    parameters.fixed.pzA = pzA; %basis choice probability (for Z basis)
    parameters.fixed.pzB = pzB;
    parameters.fixed.pxB = pxB;
    parameters.fixed.pd = pd;%1e-6; %dark count probability
    parameters.fixed.f = 1;
    parameters.fixed.fullstat = 1;
    parameters.fixed.loss = loss;%0.998;
%     parameters.fixed.etad = 1; %0.045; %detector efficiency
    parameters.fixed.decoys = decoys;%, 0];%[0.48, 0.1, 0.001]; %first is signal, other two are decoys
    parameters.fixed.ext = true;

    %%%%%%%%%%%%%%%% 3.optimizable parameters %%%%%%%%%%%%%%%%
    %optional; declaring optimizable parameters automatically invokes local search optimizers
    %must be in the format of [lowerBound, initialValue, upperBound]
    
%     parameters.optimize.mu1 = [0.1,0.2,0.3]; %can optionally optimize laser intensity (need to comment it out mu1 in 'fixed' category)
    
    
end

%set the running options for the solvers
function solverOptions=setOptions()
    %%%%%%%%%%%%%%%%% global setting %%%%%%%%%%%%%%%%%
    solverOptions.globalSetting.cvxSolver = 'mosek';
    solverOptions.globalSetting.cvxPrecision = 'high';
    
    %output level:
    %0.output nothing (except errors)
    %1.output at each main iteration
    %2.output at each solver 1 FW iteration
    %3.show all SDP solver output
    solverOptions.globalSetting.verboseLevel = 2; 
    
    %%%%%%%%%%%%%%%%% parameter optimizer setting %%%%%%%%%%%%%%%%%
    solverOptions.optimizer.name = 'coordinateDescent'; %choose between 'coordinateDescent' and 'bruteForce'
    solverOptions.optimizer.linearResolution = 3; %resolution in each dimension (for brute force search and coordinate descent)
    solverOptions.optimizer.maxIterations = 2; %max number of iterations (only for coordinate descent)
    solverOptions.optimizer.linearSearchAlgorithm = 'iterative'; %choose between built-in 'fminbnd' and custom 'iterative' algorithms for linear search (only for coordinate descent)
    solverOptions.optimizer.iterativeDepth = 2; %choose depth of iteration levels; function count = depth * linearResolution (only for coordinate descent and if using 'iterative')
    solverOptions.optimizer.maxSteps = 10; %max number of steps (only for gradient descent and ADAM)
    solverOptions.optimizer.optimizerVerboseLevel = 1; %0:only output optimized result; 1:output a progress bar; 2:output at each function call

    %%%%%%%%%%%%%%%%% step1Solver setting %%%%%%%%%%%%%%%%%
    solverOptions.solver1.name = 'asymptotic_inequality';
    
    %options mainly affecting performance
    solverOptions.solver1.maxgap = 1e-6; %1e-6 for asymptotic, 2.5e-3 for finite;
    solverOptions.solver1.maxiter = 20;%10;
    solverOptions.solver1.initmethod = 1; %minimizes norm(rho0-rho) or -lambda_min(rho), use method 1 for finite size, 2 for asymptotic v1
    
    %default options
    solverOptions.solver1.linearconstrainttolerance = 1e-10;
    solverOptions.solver1.linesearchprecision = 1e-20;
    solverOptions.solver1.linesearchminstep = 1e-3;
    solverOptions.solver1.maxgap_criteria = false; %true for testing gap, false for testing f1-f0
    solverOptions.solver1.removeLinearDependence = 'none'; %'qr' for QR decomposition, 'rref' for row echelon form, empty '' for not removing linear dependence
    

    %%%%%%%%%%%%%%%%% step2Solver setting %%%%%%%%%%%%%%%%%
    solverOptions.solver2.name = 'asymptotic_inequality';
    solverOptions.solver2.epsilon = 0;
    solverOptions.solver2.epsilonprime = 1e-12;
    
    
end