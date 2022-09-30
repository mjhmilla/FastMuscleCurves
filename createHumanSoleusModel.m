function defaultHumanSoleusModel = createHumanSoleusModel(...
                                      normPevkToActinAttachmentPoint,...
                                      normMaxActiveTitinToActinDamping,...
                                      normFiberLengthAtOneNormPassiveForce,...  
                                      ecmForceFraction,...
                                      useWLCTitinModel,...
                                      useCalibratedCurves,...
                                      useTwoSidedTitinCurves,...
                                      smallNumericallyNonZeroNumber,...
                                      flag_enableNumericallyNonZeroGradients,...
                                      scaleOptimalFiberLength,...
                                      scaleMaximumIsometricTension,...
                                      passiveForceLengthCurveSettings,...
                                      flag_useOctave)

%Potential variables to expose
rigidTendonReferenceModel             = [];
elasticTendonReferenceModel           = [];



% The human soleus model is being used to simulate the titin kinematics observed
% by Trombitas et al. during a passive stretch. As this experiment is no way
% dependent on the architectural properties of the human soleus I have not gone
% to the effort to fit the active and passive properties of the model to data.

humanSoleusActiveForceLengthData  = [];
humanSoleusPassiveForceLengthData = [];

flag_solveForOptimalFiberLengthOfBestFit         = 0; 
shiftLengthActiveForceLengthCurveDescendingCurve = 0.;

[humanSoleusMusculotendonProperties, ...
 humanSoleusSarcomereProperties,...
 humanSoleusActiveForceLengthData,...
 humanSoleusPassiveForceLengthData] = ...
    createHumanSoleusParameters( ...
                        scaleOptimalFiberLength,...
                        scaleMaximumIsometricTension,...
                        normFiberLengthAtOneNormPassiveForce,...
                        normPevkToActinAttachmentPoint,...
                        normMaxActiveTitinToActinDamping,...
                        ecmForceFraction,...
                        flag_useOctave);

createMusculoTendonFcn = ...
  @(argScaleFiberLength,argScaleFiso)createHumanSoleusParameters(...
                                        argScaleFiberLength,...
                                        argScaleFiso,...
                                        normFiberLengthAtOneNormPassiveForce,...
                                        normPevkToActinAttachmentPoint,...
                                        normMaxActiveTitinToActinDamping,...
                                        ecmForceFraction,...
                                        flag_useOctave); 
                                        

[humanSoleusNormMuscleCurvesDefault,...
 humanSoleusMusculotendonPropertiesDefault,...
 humanSoleusSarcomerePropertiesDefault,... 
 activeForceLengthCurveAnnotationPoints,...
 humanSoleusActiveForceLengthDataDefault,...
 humanSoleusPassiveForceLengthDataDefault,...
 humanSoleusPassiveForceLengthCurveSettings]= ...
    createFittedMuscleCurves( ...
      humanSoleusMusculotendonProperties,...
      humanSoleusSarcomereProperties,...
      useWLCTitinModel,...
      useCalibratedCurves,...
      useTwoSidedTitinCurves,...
      humanSoleusActiveForceLengthData,...
      humanSoleusPassiveForceLengthData,...
      passiveForceLengthCurveSettings,...
      shiftLengthActiveForceLengthCurveDescendingCurve,...
      flag_enableNumericallyNonZeroGradients,...
      smallNumericallyNonZeroNumber,...
      flag_solveForOptimalFiberLengthOfBestFit,...
      createMusculoTendonFcn,...      
      flag_useOctave);



defaultHumanSoleusModel = struct('musculotendon',...
                            humanSoleusMusculotendonPropertiesDefault,...
                            'sarcomere',...
                            humanSoleusSarcomerePropertiesDefault,...
                            'falData',...
                            humanSoleusActiveForceLengthDataDefault,...
                            'fpeData',...
                            humanSoleusPassiveForceLengthDataDefault,...
                            'curves',...
                            humanSoleusNormMuscleCurvesDefault,...
                            'fitting',...
                            []);
                      



   
