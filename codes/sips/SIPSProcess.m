function out = SIPSProcess(S1,procStruct,Q,C)
% Here we attempt to repeat the processing scripts in PS process, however
% with only a single input states

% the only two mandatory arguments
fwx = procStruct.fwx;
dz = procStruct.dz;
dopTh = 0.7;
wcorr = [];
rc = [];
dzres = 4.8;
cumulative = false; % output cumulative signal
unwrapMask = [];
noDepthCorr = true;
roiz = 1:1000; % axial range to consider
fwaxial = 1;% for comparison without spectral binning
fwz = 5; % axial filtering range of local ret
fw = 30; % lateral filtering to estimate surface signal and rotation due to sheath birefringence
dzres = 4.8;
rc = [];
wcorr = [];
cumulative = true; % flag for processing cumulative signal
sysComp = systemCompensation;
outputLevel = 2;
useSVD = 1;


fnames = fieldnames(procStruct);
for ind = 1:numel(fnames)
    if strcmp(fnames{ind},'fwz')
        fwz = procStruct.fwz;
    elseif strcmp(fnames{ind},'rc')
        rc = procStruct.rc;
    elseif strcmp(fnames{ind},'wcorr')
        wcorr = procStruct.wcorr;
    elseif strcmp(fnames{ind},'dzres')
        dzres = procStruct.dzres;
    elseif strcmp(fnames{ind},'dopTh')
        dopTh = procStruct.dopTh;
    end
end

dim = size(S1);
% manage case of no spectral binning
if numel(dim)<4 % no spectral binning used
    S1 = permute(S1,[1,2,4,3]);% introduce third dimension
    dim = size(S1);
end

if dim(4) == 3 % 3-component Stokes vector was provided
    S1 = cat(4,sqrt(dot(S1,S1,4)),S1);
end

nx = (round(fwx*1.5)-1)/2;
nx = linspace(-nx,nx,round(fwx*1.5))*2*sqrt(log(2))/fwx;
h = exp(-nx.^2);
if fwz>1
    nz = (round(fwaxial*1.5)-1)/2;
    nz = linspace(-nz,nz,round(fwaxial*1.5))*2*sqrt(log(2))/fwaxial;
    h = exp(-nz(:).^2)*h;
end
h = h/sum(h(:));
S1f = imfilter(S1,h,'circular');

% final normalization
QUVf = sqrt(dot(S1f(:,:,:,2:4),S1f(:,:,:,2:4),4));
dop = mean(QUVf./S1f(:,:,:,1),3);

% Make Stokes into Jones vector for GDA
tom = Stokes2Jones(permute(S1f(:,:,:,2:4),[4,1,2,3]));

redFact = 10e4;
%tom = tom ./ redFact;
tom = tom./sqrt(abs(tom(1,:,:,:)).^2+abs(tom(2,:,:,:)).^2);
Nbins = size(S1,3);

% Estimate C & Q
%[C,Q] = estimateCQSIPS(S1f);

dimQ = size(Q);

if dimQ(1) ~= 2
    %Q = reshape(Q,[2,2,Nbins]);
    C = reshape(C,[2,2,Nbins]);

end

% Dimensions & parameters

% dim = 150;
% startRow = 1;
% startCol = 250;
% 
% tom = tom(:,startRow:startRow+dim,startCol:startCol+dim,:);
dimA = size(tom,2);
dimL = size(tom,3);

maxIter = 60;
tol = 1e-3;
alpha = 0.8;
gamma = 0.2;
step = 0.4;


if useSVD == 0 
    [outTom,gdTime] = vectorizedMomentumGDA(C,Q,tom,dimA,dimL,maxIter,tol,alpha,gamma,dop);
    r = cat(1,outTom,zeros(1,dimA,dimL));
else
    [r,singVals]= svdDecomposition(C,Q,tom);
end


MM = squeeze(makeRot(r));
MM = euclideanRotation(MM + MM([1,4,7,2,5,8,3,6,9]',:,:).*[1;1;-1;1;1;-1;-1;-1;1],true);
%MM = reshape(MM,[9,size(MM,3),size(MM,4)]);
%size(MM)
dmn = MatrixMultiply(MM(:,2:end,:),MM([1,4,7,2,5,8,3,6,9],1:end-1,:));
locoa = cat(2,zeros(3,1,size(MM,3)),decomposeRot(dmn));

Omegaf = imfilter(permute(locoa,[2,3,1]),ones(dz,1)/dz);
retFinal = sqrt(sum(Omegaf.^2,3))/dzres/pi*180*100;

% Average
% for indw = 1:9
%     r_bin = cat(1,tom(:,:,:,indw),zeros(1,dimA,dimL));
%     MMav(:,:,:,indw) = makeRot(r_bin) ;
% end
% MMav = mean(MMav,4);
% 
% dmn_av = MatrixMultiply(MMav(:,2:end,:),MMav([1,4,7,2,5,8,3,6,9],1:end-1,:));
% locoa_av = cat(2,zeros(3,1,size(MMav,3)),decomposeRot(dmn_av));
% 
% Omegaf_av = imfilter(permute(locoa_av,[2,3,1]),ones(dz,1)/dz);
% retAvFinal = sqrt(sum(Omegaf_av.^2,3))/dzres/pi*180*100;

out = struct;
out.dop = sum(dop,3);
out.tom = tom;
out.MM = MM;
out.omega = Omegaf;
out.locoa = locoa;
out.ret = retFinal;
out.S1f = S1f;
out.singVals = singVals;

end