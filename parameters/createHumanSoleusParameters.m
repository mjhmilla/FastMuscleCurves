function [humanSoleusMusculotendonProperties, ...
          humanSoleusSarcomereProperties,...
          humanSoleusActiveForceLengthData,...
          humanSoleusPassiveForceLengthData] = ...
          createHumanSoleusParameters(...
            scaleOptimalFiberLength,...
            scaleMaximumIsometricTension,...
            normFiberLengthAtOneNormPassiveForce,...
            normPevkToActinAttachmentPoint,...
            normMaxActiveTitinToActinDamping,...
            ecmForceFraction,...
            flag_useOctave)
%%
% This function uses data from the literature to return a series of structs that
% contain the necessary musculotendon properties, sarcomere properties, and
% processed data to construct Opus 31 for a human soleus muscle which is a 
% rather detailed muscle model.
% Please review the script for details and references to the literature.
% Default values that have been determined by simulation have been noted
% as such in the script
%
% @param scaleOptimalFiberLength: scales the optimal fiber length
%
% @param scaleMaximumIsometricTension: scales the maximum isometric tension
%
% @param flag_useOctave
%   Setting this to 1 will ensure that no parts of the code that are 
%   incompatible with octave are called.
%
% @return Four structs:
%   humanSoleusMusculotendonProperties
%     architectural properties and some gross mechanical properites
%   humanSoleusSarcomereProperties
%     lengths of all of the various filaments of a human soleus along with 
%     detailed information regarding the segment lengths of titin
%   humanSoleusActiveForceLengthData
%     normalized versions of data used to to fit the active-force-length 
%     curve (presently empty).
%   humanSoleusPassiveForceLengthData
%     normalized versions of data used to to fit the passive-force-length 
%     curve (presently empty).
%%
%scaleOptimalFiberLength               = 1.0; %user-settable parameter
%scaleMaximumIsometricTension          = 1.0; %user-settable parameter

%%
% This is the normalized fiber length at which the extrapolated 
% passive-force-length curve is expected to develop 1 normalized force.
% Here this is the passive force length curve that is fitted to the data
% of Herzog & Leonard 2002.
%%

%Get the default sarcomere properties for a human soles                               
[humanSoleusSarcomereProperties] =...
  getMammalianSkeletalMuscleNormalizedSarcomereProperties(...
    'human',...
    normFiberLengthAtOneNormPassiveForce,...
    normPevkToActinAttachmentPoint,...
    normMaxActiveTitinToActinDamping,...
    ecmForceFraction,...
    []);

%%
%  As is typical when simulating young adults we will use 10 lopt/s as the
%  maximum velocity. It is entirely possible that this value is high because
%  Hill models are typically used in simulation. During movements that 
%  involve high eccentric loading (like running) a Hill model will look
%  weaker than the proposed model. A way to compensate for this is to
%  make the muscle faster. 
%
%  Hopefully at some point in-vivo fiber kinematics will be measured during
%  fast contractions and this number can be updated.
%
%  Thelen DG. Adjustment of muscle mechanics model parameters to simulate 
%  dynamic contractions in older adults. J Biomech Eng 2003;125(1):70–7.
%
%  Arnold EM, Delp SL. Fibre operating lengths of human lower limb muscles
%  during walking. Philos Trans R Soc Lond B, Biol Sci 2011;366(1570):1530–9.
%%

% Update: 22 April 2024
%
% When a Hill force-velocity curve is fit to the data of Fontana et al.
% vceMax is 5.66 l/s. The minimum fiber velocity that is consistent
% with Nikolaidou et al.'s data is 5.57 l/s. The average is 5.615 lo/s.
%
% The value of the force velocity curve at half vmax is 0.21873 and 0.1975
% respectively, the average of which is 0.2081.
%
% de Brito Fontana H, Roesler H, Herzog W. In vivo vastus lateralis 
% force–velocity relationship at the fascicle and muscle tendon unit level. 
% Journal of Electromyography and Kinesiology. 2014 Dec 1;24(6):934-40.
%
% Nikolaidou ME, Marzilger R, Bohm S, Mersmann F, Arampatzis A. Operating 
% length and velocity of human M. vastus lateralis fascicles during 
% vertical jumping. Royal Society Open Science. 2017 May 3;4(5):170185.
%

maximumNormalizedFiberVelocity = 5.615; % in units of norm fiber lengths/second
forceVelocityMultiplierAtHalfMaximumFiberVelocity = 0.2081;  

% Update: 22 April 2024

% The slow-twitch fibers plotted in Fig. 3 of Ranatunga 1984 develop a 
% normalized force of 0.1 at half the maximum contraction velocity. I am
% assuming that the rats slow-twich normalized force-velocity curve will also
% be similar for human, as they are bot mammals. [To do: look for better numbers]
%
% Ranatunga KW. The force‐velocity relation of rat fast‐and slow‐twitch muscles 
% examined at different temperatures. The Journal of physiology. 1984 Jun 1;
% 351(1):517-29.
% forceVelocityMultiplierAtHalfMaximumFiberVelocity = 0.1;  


% Fit to in-vivo data of the human Achilles tendon from Magnusson et al.
tendonStrainAtOneNormForce      = 0.049; 

%Get the (formatted) experimental data on the active/passive
%force-length curves
useElasticTendonExp = 1;


normPlateauOffset = ...
  humanSoleusSarcomereProperties.normMyosinBareHalfLength;


%Get the default musculotendon properties for the human soleus
[humanSoleusMusculotendonProperties] = ...
  getHumanSoleusMusculotendonProperties(...
            maximumNormalizedFiberVelocity,...
            forceVelocityMultiplierAtHalfMaximumFiberVelocity,...
            tendonStrainAtOneNormForce,...
            scaleOptimalFiberLength,...                              
            scaleMaximumIsometricTension,...
            normPlateauOffset,...
            useElasticTendonExp,...
            flag_useOctave);
          

humanSoleusActiveForceLengthData = [];
humanSoleusPassiveForceLengthData= [];

  

