function out = SIPSOAProcess(S1,procStruct,Q,C)
% Here we attempt to repeat the processing scripts in PS process, however
% with only a single input states

% the only two mandatory arguments
fwx = procStruct.fwx;
dz = procStruct.dz;
fwz = 1;
dopTh = 0.6;
wcorr = [];
rc = [];
dzres = 4.8;
outputLevel = 1;
cumulative = true; % output cumulative signal
unwrapMask = [];
noDepthCorr = true;
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

%tic
if useSVD == 0 
    [outTom,gdTime] = vectorizedMomentumGDA(C,Q,tom,dimA,dimL,maxIter,tol,alpha,gamma,dop);
    r = cat(1,outTom,zeros(1,dimA,dimL));
    out.outTom = outTom;
    out.gdTime = gdTime;

else
    r= svdDecomposition(C,Q,tom);
end
%outTom = vectorizeGDAOrthoConstraint(C,Q,tom,dimA,dimL,maxIter,tol,alpha);
%outTom = vectorizedGDARiemannian(C,Q,tom,dimA,dimL,maxIter,tol,step,dop);

%toc

MM = squeeze(makeRot(r));
MM = euclideanRotation(MM + MM([1,4,7,2,5,8,3,6,9]',:,:).*[1;1;-1;1;1;-1;-1;-1;1],true);
dmn = MatrixMultiply(MM(:,2:end,:),MM([1,4,7,2,5,8,3,6,9],1:end-1,:));
locoa = cat(2,zeros(3,1,size(MM,3)),decomposeRot(dmn));


Omegaf = imfilter(permute(locoa,[2,3,1]),ones(dz,1)/dz);
retFinal = sqrt(sum(Omegaf.^2,3))/dzres/pi*180*100;

% Get OA Information
[sqMMinv,OAunwrapped,unwrapMask] = unwrapOAx(MM,dop,dopTh,unwrapMask);

sqMMinv = reshape(sqMMinv,[9,size(MM,2),size(MM,3)]);

dmn = MatrixMultiply(sqMMinv(:,1:end-1,:),MatrixMultiply(MM(:,2:end,:),sqMMinv(:,1:end-1,:)));
locoa = decomposeRot(dmn);
VV = MatrixMultiply(makeRot(locoa/2),MatrixMultiply(sqMMinv([1,4,7,2,5,8,3,6,9],1:end-1,:),sqMMinv(:,2:end,:)));
vv = decomposeRot(VV);

% first, use DOP to mask meaningless areas
mask = dop>dopTh;
mask = mask(1:end-1,:)&mask(2:end,:);

% clean this mask up with an additional criterion on the uniformity of
% the local optic axis over an axial range.
maskoa = abs(imfilter(exp(1i*squeeze(atan2(locoa(1,:,:),locoa(2,:,:)))),ones(fwz,1)/fwz))>.9;
 
mask = mask&maskoa;

vvphi = squeeze(vv(3,:,:));
vvphi(~mask) = 0;    
vvcorr = shiftdim(cumsum(cat(1,zeros(1,size(vvphi,2)),vvphi(1:end-1,:))),-1);

oa(1,:,:) = locoa(1,:,:).*cos(vvcorr) - locoa(2,:,:).*sin(vvcorr);
oa(2,:,:) = locoa(1,:,:).*sin(vvcorr) + locoa(2,:,:).*cos(vvcorr);
oa(3,:,:) = zeros(size(locoa(1,:,:)));

locoa = cat(2,zeros(3,1,size(locoa,3)),locoa);
if noDepthCorr% if set to true, ignore the depth-corrected signal
    oa = locoa;
else
    oa = cat(2,zeros(3,1,size(oa,3)),oa);    
end

if outputLevel>0
    if ~cumulative % without depth correction
        Omegaf = imfilter(permute(locoa,[2,3,1]),ones(fwz,1)/fwz);
        retFinal2 = sqrt(sum(Omegaf.^2,3))/dzres/pi*180*100;
        phi2 = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));
        out.phi2 = phi2;
        out.ret2 = retFinal2;
    else % outputting the cumulative signal
        W = decomposeRot(MM);
        Omegaf = permute(W,[2,3,1]);
        retFinal2 = sqrt(sum(Omegaf.^2,3))/pi*100;% scale pi to 100
        phi2 = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));
        out.phi2 = phi2;
        out.ret2 = retFinal2;
    end
    out.OAunwrapped = OAunwrapped;
end


Omegaf = imfilter(permute(locoa,[2,3,1]),ones(dz,1)/dz);
phi = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));

out.dop = dop;
out.mask = mask;
out.ret = retFinal;
out.phi = phi;
out.fwz = fwz;
out.fwx = fwx;
out.unwrapMask = unwrapMask;
out.tom = tom;
out.MM = MM;
out.omega = Omegaf;
out.locoa = locoa;
out.S1f = S1f;


