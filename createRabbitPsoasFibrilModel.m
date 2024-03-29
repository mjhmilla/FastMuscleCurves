function defaultRabbitPsoasFibril = createRabbitPsoasFibrilModel(...
                                      normCrossBridgeStiffness,...
                                      normCrossBridgeDamping,...
                                      normPevkToActinAttachmentPoint,...
                                      normMaxActiveTitinToActinDamping,...
                                      normFiberLengthAtOneNormPassiveForce,...
                                      ecmForceFraction,...
                                      titinMolecularWeightInkD,...
                                      useWLCTitinModel,...                                      
                                      useCalibratedCurves,...
                                      useTwoSidedTitinCurves,...
                                      smallNumericallyNonZeroNumber,...
                                      flag_enableNumericallyNonZeroGradients,...
                                      scaleOptimalFiberLength,...
                                      scaleMaximumIsometricTension,...
                                      flag_useOctave)

% The human soleus model is being used to simulate the titin kinematics observed
% by Trombitas et al. during a passive stretch. As this experiment is no way
% dependent on the architectural properties of the human soleus I have not gone
% to the effort to fit the active and passive properties of the model to data.
rabbitPsoasFibrilActiveForceLengthData  = [];
rabbitPsoasFibrilPassiveForceLengthData = [];


% Check the arguments
if( smallNumericallyNonZeroNumber <= 0 &&...
    flag_enableNumericallyNonZeroGradients  )
  disp('Warning: flag_enableNumericallyNonZeroGradients is set to 1 and');
  disp('         smallNumericallyNonZeroNumber <= 0');
  disp('Setting smallNumericallyNonZeroNumber to be sqrt(sqrt(eps))');
  smallNumericallyNonZeroNumber = sqrt(sqrt(eps));
end

if(isempty(normPevkToActinAttachmentPoint))
  normPevkToActinAttachmentPoint        = 0.5;

  disp('normPevkToActinAttachmentPoint empty using default:');  
  fprintf('\t%f\n',normPevkToActinAttachmentPoint);
end

if(isempty(normFiberLengthAtOneNormPassiveForce))
  normFiberLengthAtOneNormPassiveForce        = 1.367732948060934e+00;

  disp('normFiberLengthAtOneNormPassiveForce empty using default');
  fprintf('\t%f\n',normFiberLengthAtOneNormPassiveForce);
end



[rabbitPsoasFibrilMusculotendonProperties, ...
 rabbitPsoasFibrilSarcomereProperties] = ...
    createRabbitPsoasFibrilParameters(  ...
                        scaleOptimalFiberLength,...
                        scaleMaximumIsometricTension,...
                        normFiberLengthAtOneNormPassiveForce,...
                        normPevkToActinAttachmentPoint,...
                        normMaxActiveTitinToActinDamping,...
                        ecmForceFraction,...
                        titinMolecularWeightInkD,...
                        flag_useOctave);

rabbitPsoasFibrilSarcomereProperties.normCrossBridgeStiffness ...
                                    =normCrossBridgeStiffness;

rabbitPsoasFibrilSarcomereProperties.normCrossBridgeDamping ...
                                    =normCrossBridgeDamping;

createMusculoTendonFcn = ...
  @(argScaleFiberLength,argScaleFiso)createRabbitPsoasFibrilParameters(...
                                        argScaleFiberLength,...
                                        argScaleFiso,...
                                        normFiberLengthAtOneNormPassiveForce,...
                                        normPevkToActinAttachmentPoint,...
                                        normMaxActiveTitinToActinDamping,...
                                        ecmForceFraction,...
                                        flag_useOctave); 
                                        


%We have no data to fit to, and so these options cannot be used
flag_solveForOptimalFiberLengthOfBestFit         = 0; 
shiftLengthActiveForceLengthCurveDescendingCurve = 0.;

passiveForceLengthCurveSettings = [];

[rabbitPsoasFibrilNormMuscleCurvesDefault,...
 rabbitPsoasFibrilMusculotendonPropertiesDefault,...
 rabbitPsoasFibrilSarcomerePropertiesDefault,... 
 activeForceLengthCurveAnnotationPoints,...
 rabbitPsoasFibrilActiveForceLengthDataDefault,...
 rabbitPsoasFibrilPassiveForceLengthDataDefault,...
 rabbitPsoasFibrilPassiveForceLengthCurveSettings]= ...
    createFittedMuscleCurves( ...
      rabbitPsoasFibrilMusculotendonProperties,...
      rabbitPsoasFibrilSarcomereProperties,...
      useWLCTitinModel,...      
      useCalibratedCurves,...
      useTwoSidedTitinCurves,...
      rabbitPsoasFibrilActiveForceLengthData,...
      rabbitPsoasFibrilPassiveForceLengthData,...
      passiveForceLengthCurveSettings,...
      shiftLengthActiveForceLengthCurveDescendingCurve,...
      flag_enableNumericallyNonZeroGradients,...
      smallNumericallyNonZeroNumber,...
      flag_solveForOptimalFiberLengthOfBestFit,...
      createMusculoTendonFcn,...
      flag_useOctave);



defaultRabbitPsoasFibril = struct('musculotendon',...
                            rabbitPsoasFibrilMusculotendonPropertiesDefault,...
                            'sarcomere',...
                            rabbitPsoasFibrilSarcomerePropertiesDefault,...
                            'falData',...
                            rabbitPsoasFibrilActiveForceLengthDataDefault,...
                            'fpeData',...
                            rabbitPsoasFibrilPassiveForceLengthDataDefault,...
                            'curves',...
                            rabbitPsoasFibrilNormMuscleCurvesDefault,...
                            'fitting',...
                            []);
              


   
