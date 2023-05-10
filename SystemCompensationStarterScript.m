%% Generating your own system compensation file
% Requires: Significant heterogenous data (>50 B-scans)
%           Polarization Diverse Detection

% First, obtain 50+ B-scans of spectrally-binned Stokes vectors from your system
% (feel free to use our MAKESTOKES(TOM,DIM))
% where DIM refers to a 3 or 4 component Stokes vector. We will be using
% 3-component stokes vectors

S1 = %... add your own processing scheme

% Place Stokes vectors in the following structure format (5D single)
% S1 = [axial dimension, lateral dimension, # spectral bins, 3 (stokes parameters), # slices]


% Define processing parameters
syscomstruct.isMod = 0; % is this data using A-line polarization modulation 
syscomstruct.isFilt = 0; % is this data filtered prior


% Feed into processing pipeline
sysComp = estimateSystemCompensation(S1,syscomstruct);

% Evaulate System Compensation
[ill_dgd, det_dgd,sps] = characterizeSysCom(sysComp);
%%
disp("Spectral Polarization Spread = "+sps)
disp("")
disp("         If SPS < 1, single input processing may have higher error (see [1])")
disp("         You can consider adding extra PMD through a small piece of PM fiber or adding a polarization controler")

