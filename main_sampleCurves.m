clc;
close all;
clear all;

addpath(fullfile('curves','BezierCurveLibrary'));

%%
% In output/tables/curves you'll find these 4 directories that contain
% Bezier coefficients for muscle model curves:
%
% QuadraticBezierHumanCurves
% QuadraticBezierFelineCurves
% QuadraticBezierCurvesOne
% QuadraticBezierCurvesZero
%
%
%
% The human and feline curves are in the directories that you would expect.
% It is worth noting that the 
%
% ...forceLengthProximalTitinCurve
% ...forceLengthDistalTitinCurve
%
% have been formed assuming that one specific location of the PEVK segment
% binds to actin. This location is set in
%
%   main_createExplicitBezierSplineMuscleCurves.m 
%
% by the variable 
%
%   normPevkToActinAttachmentPointDefault on line 103 
%
% (currently set to 0.5) where a value of 0 corresponds to the most
% proximal part of the PEVK segment binding to actin (which makes titin
% compliant during active lengthening) and a value of 1 corresponds to the
% most distal part of the PEVK segment binding to actin (which makes titin
% quite stiff). 
% 
% If you want to be cautious, then start with values of 
% 0 for normPevkToActinAttachmentPointDefault and only increase if it is
% needed to fit experimental data. Here a value of 0.5 has been used to 
% allow the cat soleus model to fit Herzog and Leonard 2002. Is this
% appriopriate for your simulations? I don't know. This is cutting edge
% research at the moment, so you'll have to do your own detective work.
%
% Herzog W, Leonard TR. Force enhancement following stretching of skeletal 
% muscle: a new mechanism. Journal of Experimental Biology. 
% 2002 May 1;205(9):1275-83.
%%

%%
% If you want to sample curves from a different directory, then update the
% 'QuadraticBezierHumanCurves' and 'QuadraticBezierHumanCurvesSampled'
% below
%%
curveFolder = fullfile('output','tables','curves',...
                    'QuadraticBezierHumanCurves');
sampleFolder = fullfile('output','tables','curves',...
                    'QuadraticBezierHumanCurvesSampled');

%Reads in the Bezier coefficients
structOfMuscleCurves = readBezierCurveStructuresFromCSV([curveFolder,filesep]);
nPts                = 100;
muscleName          = '';
flag_extendedDomain = 0; 

status = writeQuadraticBezierCurveSamplesToCSV(...
            structOfMuscleCurves,nPts,muscleName,[sampleFolder,filesep],...
            flag_extendedDomain);







