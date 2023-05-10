%% Single Input Polarization Sensitive (SIPS) Processing of OFDI Data
% This is a demonstration of the SIPS processing pipeline of non-modulating
% OFDI data acquired with polarization diverse detection

% First, we must generate the Stokes vectors for our measurements
% This can be done by obtaining the processed, complex, spectrally-binned [2] tomogram from your system
% and running MAKESTOKES(TOM,DIM)
% 
% [1] G.L. Jones, Q. Xiong, X. Liu, B. E. Bouma, and M. Villiger. "Single-Input Polarization-Sensitive Optical Coherence Tomography Through
% a Catheter", arxiv.com

% Load Tomogram
load('examples\savedTom1.mat')

% Load System Compensation
% To make your own system compensation, see SystemCompensationStarterScript.m
load('examples\SystemCompensation.mat') 

% Chracterize PMD
[ill_dgd, det_dgd,sps] = characterizeSysCom(sysComp);
% if sps < 1.0, potentialy high overall error of measurments (see [1]) 

% Find Symmetric and Asymmetric Corrctions
N = size(sysComp.alignRotVec,2); % # spectral bins
Q = reshape(makeJones(-sysComp.alignRotVec),[4,N]); % Symmetric Compensation
C = reshape(makeJones(sysComp.symRotVec),[4,N]); % Asymmetric Compensation
dopThresh = 0.7;

% Number of Cross Sections to Process
n_slices = size(binned_tom,2);


%% Process Sections
% Processing structure
pstruct.fwx = 6; % lateral filtering
pstruct.dz = 5; % axial filtering
pstruct.dzres = 4.8; % axial resolution
for slice_ind = 1:n_slices

    % Use full spectrum tomogram to create intensity image
    tomogram = full_tom{slice_ind};
    intensity = abs(squeeze(tomogram(1,:,:))).^2 + abs(squeeze(tomogram(2,:,:))).^2;

    % Observe intensity image
    figure,
    subplot(2,2,1)
    imagesc(10*log10(intensity)),colormap('gray')
    title("Intensity Image (dB)")

    % Convert binned vectors to stokes for processing
    binned_stokes = makeStokes(binned_tom{slice_ind},3);

    % Send into SIPS Processing Pipeline
    % Retardance & DOP ONLY
    outRet = SIPSProcess(binned_stokes,pstruct,Q,C);

    % Retardance, DOP & Optic Axis
    out = SIPSOAProcessCath(binned_stokes,intensity,pstruct,Q,C);
    
    % Plot Images
    dopImage = outRet.dop;
    retImage = outRet.ret;
    retImage(dopImage<dopThresh) = 0;
    oaImage = out.phi;
    oaImage(dopImage<dopThresh) = 0;

    subplot(2,2,2)
    imagesc(retImage, [0,120])
    title("Retardance")
    subplot(2,2,3)
    imagesc(outRet.dop,[0,1])
    title("Degree of Polarization")
    subplot(2,2,4)
    imagesc(oaImage,[-pi,pi])
    title("Depth-Resolved OA")
    pause
end