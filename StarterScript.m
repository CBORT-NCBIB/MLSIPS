%% Single Input Polarization Sensitive (SIPS) Processing of OFDI Data
% This is a demonstration of the SIPS processing pipeline of non-modulating
% OFDI data acquired with polarization diverse detection

% First, we must generate the Stokes vectors for our measurements
% This can be done by obtaining the processed, complex, spectrally-binned [2] tomogram from your system
% and running MAKESTOKES(TOM,DIM)
% 
% [1] G.L. Jones, Q. Xiong, X. Liu, B. E. Bouma, and M. Villiger. "Single-Input Polarization-Sensitive Optical Coherence Tomography Through
% a Catheter", arxiv.com

%% Add supporting code and example data to path
addpath(genpath(fullfile('codes')));
addpath(genpath(fullfile('examples')));

%% Load Tomogram
% Loaded inputs are as follows
% full_tom: 
%   1x1 cell containing a 2 x 1024 x 512 complex single value matrix. This
%   contains the unbinned complex reconstructed tomogram
%   2 - tomogram channels, 1024 - 1D depth scan size, 512 - # of 1D scans per lateral image
% binned_tom: 
%   1x1 cell containing a 2 x 1024 x 512 x 9 complex single value matrix. This
%   contains the binned complex reconstructed tomogram. 
%   2 - tomogram channels, 1024 - 1D depth scan size, 512 - # of 1D scans
%   per lateral image, 9 - number of bins
load('examples\savedTom.mat')


%% Load System Compensation
% To make your own system compensation, see SystemCompensationStarterScript.m
% Loaded input includes 1 systemCompensation structure. 
% Notable parameters include
    % symRotVec -> rotation vector which describes asymmetric compensation
    % alignRotVec -> rotation vector which describes symmetric compensation
    % visualize -> shows value of the above across spectral bins
load('examples\SystemCompensation.mat') 

%% Load Colormaps
% Obtain 3 colormaps which help visualize the depolarization (D),
% birefringence(B) and optic axis (OA)
load('examples\colormaps.mat') 

%% Chracterize PMD
sps = getSPS(sysComp);
% if sps < 1.0, potentialy high overall error of measurments (see [1]) 

%% Find Symmetric and Asymmetric Corrctions
N = size(sysComp.alignRotVec,2); % Number of spectral bins
Q = reshape(makeJones(-sysComp.alignRotVec),[4,N]); % Symmetric Compensation Jones Matrix
C = reshape(makeJones(sysComp.symRotVec),[4,N]); % Asymmetric Compensation Jones Matrix
dopThresh = 0.7;

%% Number of Cross Sections to Process
n_slices = size(binned_tom,2);


%% Process Sections
% Processing structure
pstruct.fwx = 6; % lateral filtering (in px)
pstruct.dz = 5; % axial filtering (in px)
pstruct.dzres = 4.8; % axial resolution (in um)

for slice_ind = 1:n_slices

    % Use full spectrum tomogram to create intensity image
    tomogram = full_tom{slice_ind};
    intensity = abs(squeeze(tomogram(1,:,:))).^2 + abs(squeeze(tomogram(2,:,:))).^2;

    % Convert binned vectors to stokes for processing
    disp("Converting to Stokes...")
    binned_stokes = makeStokes(binned_tom{slice_ind},3);

    % Send into SIPS Processing Pipeline
    % Birefringence & Depolarization ONLY
    disp("Computing Birefringence & Depolarization...")
    outRet = SIPSProcess(binned_stokes,pstruct,Q,C);

    % Birefringence & Depolarization & Optic Axis
    disp("Computing Birefringence, Depolarization & Optic Axis...")
    out = SIPSOAProcessCath(binned_stokes,intensity,pstruct,Q,C);
    
    %% Plot Images
    dopImage = outRet.dop;
    retImage = outRet.ret;
    retImage(dopImage<dopThresh) = 0;
    oaImage = out.phi;
    oaImage(dopImage<dopThresh) = 0;


    figure,
    ax1= subplot(2,2,1)
    imagesc(10*log10(intensity)),
    title("Intensity Image (dB)")
    colormap(ax1,'gray')
    ax2 = subplot(2,2,2)
    imagesc(retImage, [0,100])
    title("Birefringence")
    colormap(ax2,cmapB)
    ax3 = subplot(2,2,3)
    imagesc(1-outRet.dop,[0,0.4])
    title("Depolarization")
    colormap(ax3,cmapD)
    ax4 = subplot(2,2,4)
    imagesc(oaImage,[-pi,pi])
    title("Depth-Resolved OA")
    colormap(ax4,cmapOA)
end