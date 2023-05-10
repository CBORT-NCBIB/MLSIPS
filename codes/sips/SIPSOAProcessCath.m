function out = SIPSOAProcessCath(S1,int,procStruct,Q,C)
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

tic

if useSVD == 0 
    [outTom,gdTime] = vectorizedMomentumGDA(C,Q,tom,dimA,dimL,maxIter,tol,alpha,gamma,dop);
    r = cat(1,outTom,zeros(1,dimA,dimL));
else
    if outputLevel>1
    [r,singVals]= svdDecomposition(C,Q,tom);
    out.singVals = singVals;
    else
    r= svdDecomposition(C,Q,tom);
    end
end
%outTom = vectorizeGDAOrthoConstraint(C,Q,tom,dimA,dimL,maxIter,tol,alpha);
%outTom = vectorizedGDARiemannian(C,Q,tom,dimA,dimL,maxIter,tol,step,dop);

toc

MM = squeeze(makeRot(r));

if outputLevel > 1
    out.OA1 = decomposeRot(MM);
end

%MM = euclideanRotation(MM + MM([1,4,7,2,5,8,3,6,9]',:,:).*[1;1;-1;1;1;-1;-1;-1;1],true);
dmn = MatrixMultiply(MM(:,2:end,:),MM([1,4,7,2,5,8,3,6,9],1:end-1,:));
locoa = cat(2,zeros(3,1,size(MM,3)),decomposeRot(dmn));
Omegaf = imfilter(permute(locoa,[2,3,1]),ones(dz,1)/dz);
retFinal = sqrt(sum(Omegaf.^2,3))/dzres/pi*180*100;
%% detect catheter
cath = findCatheter(int);

if outputLevel > 0
    out.cath = cath;
end

% identify signal from catheter
%Itotf = I1f + I2f;
% Take intensity weighted average between first and second peak
mm = (bsxfun(@minus,(1:size(dop,1))',cath(1,:)-10)>0)&(bsxfun(@minus,(1:size(dop,1))',cath(2,:)+10)<0);
rrz = max(min(cath(1,:))-10,1):max(cath(2,:))+10;
mm = mm(rrz,:);
ww = shiftdim(repmat(max(dop(rrz,:)*2-1,0).*mm,[1,1,size(MM,4)]),-1);
Mm = squeeze(bsxfun(@rdivide,sum(bsxfun(@times,MM(:,rrz,:,:),ww),2),sum(ww*9,2)));
Mm = imfilter(Mm,ones(1,fw)/fw,'circular');
Mm = euclideanRotation(Mm);
WFit = decomposeRot(Mm);

if outputLevel > 0
    out.surf = WFit;
end

%% unwrap fit along the B-scan
ss = repmat(cat(2,0,mod(cumsum(sqrt(sum(diff(WFit,[],2).^2,1))>pi,2),2)),[3,1])>0;
delta = 2*pi*bsxfun(@rdivide,WFit,sqrt(sum(WFit.^2,1)));

WFitc = WFit;
WFitc(ss) = WFit(ss) - delta(ss);

% Check if overall a subtraction/addition of pi would make more sense
n = -1:1;
for ind = 1:3
    WW = WFitc + 2*pi*n(ind)*bsxfun(@rdivide,WFitc,sqrt(sum(WFitc.^2,1)));
    meanRet(ind) = mean(mean(sqrt(sum(WW.^2,1))));
end

[~,mp] = min(meanRet);
WW = WFitc + 2*pi*n(mp)*bsxfun(@rdivide,WFitc,sqrt(sum(WFitc.^2,1)));
WW(3,:,:) = 0;

if outputLevel > 0
    out.surfCorr = WW;
end

ball = decomposeBallLens(WW);% false triggers the correction of only the linear rotation part
MCorr = permute(ball.Mcorr,[1,3,2]);
MCorrT = bsxfun(@times,MCorr([1,4,7,2,5,8,3,6,9],:,:),[1;1;-1;1;1;-1;-1;-1;1]);
MM = MatrixMultiply(repmat(MCorrT,[1,size(MM,2),1]),MatrixMultiply(MM,repmat(MCorr,[1,size(MM,2),1])));

if outputLevel > 0
    out.ball = ball;
end
if outputLevel > 1
    out.OA3 = decomposeRot(MM);
end

%% now, extract signal of the outher sheath interface
% rrz = min(cath(3,:))-9:max(cath(3,:)) + 9;
% 
% ww = shiftdim(dop(rrz,:),-1);
% Mm = squeeze(bsxfun(@rdivide,sum(bsxfun(@times,MM(:,rrz,:),ww),2),sum(ww*9,2)));
% thcorr = decomposeRot(euclideanRotation(imfilter(Mm,ones(1,fw)/fw,'circular')));
% sheathRot = atan2(thcorr(2,:),thcorr(1,:));
% sheathRet = sqrt(sum(thcorr.^2,1));
% 
% if outputLevel > 0
%     out.sheath = sheathRot(:);
%     out.sheathRet = sheathRet(:);
% end
% 

% extract surface signal
mask = dop>0.8;
mask(1,:) = false; % make sure, at least the very first line is zero
%mask(1:refInds(3),:) = false; % make sure, at least the very first line is zero
mask(end,:) = false; % as well as the last one, to have an even pair of up and downward edges in each A-line

flipmask = flip(mask,1);
temp = cumsum(flipmask,1);
dm = cat(1,zeros(1,size(mask,2)),diff(flipmask,[],1));
lInds1 = find(dm>0);
lInds2 = find(dm<0);
ss = zeros(size(temp));
ss(lInds2) = temp(lInds2)-temp(lInds1);
[~,mp] = max(flip(ss,1),[],1);

% improve the surface detection by medfiltering and a 'smart' approach to
% pick the original and filtered surface signal
mpf = medfilt1(cat(2,mp,mp),fw);
stdmp = sqrt(medfilt1((cat(2,mp,mp)-mpf).^2,fw));
stdmp = circshift(stdmp(dim(2)/2 + (1:dim(2))),dim(2)/2);
mpf = circshift(mpf(dim(2)/2 + (1:dim(2))),dim(2)/2);
mp(abs(mp-mpf)>2*stdmp) = mpf(abs(mp-mpf)>2*stdmp);
mp = round(mp);
mp = max(mp + 8,cath(3,:));
%mp(mp<refInds(3)) = refInds(3);

inds = sub2ind(size(MM),ones(1,numel(mp)),mp,1:size(MM,3));
inds = bsxfun(@plus,inds,(0:size(MM,1)-1)');
Mm = MM(inds);
OAtissue = decomposeRot(euclideanRotation(imfilter(Mm,ones(1,fw)/fw,'circular')));

out.tissueAngle = atan2(OAtissue(2,:),OAtissue(1,:));
out.tissueRet = sqrt(sum(OAtissue.^2,1));
out.tissueSurf = mp;

[ang,err] = estimateOrientation(OAtissue,1);% mean 
out.sheathAngle = ang;
out.errSheathAngle = err;

% % use the signal of the sample surface to determine the best V-rotation
% % correction strategy
% [Vcorr,OAtissuep,ang,err,psi,type,ballErr,OAmodel] = refineVcorrection(ball,OAtissue);

%out.sheathAngle = ang;
%out.errSheathAngle = err;
%out.refineType = type;
%out.tissueOffset = psi;
%out.Vcorr = Vcorr;
%out.tissueAngle = squeeze(atan2(OAtissuep(2,:,:),OAtissuep(1,:,:)));
%out.tissueRet = squeeze(sqrt(sum(OAtissuep.^2,1)));
%out.ballErr = ballErr;
%out.OAmodel = OAmodel;

%Vcorr = Vcorr(2,:);% let's use the trace lsq solution
%%
% generate mask to exclude points with low DOP from correcting the local
% retardation
mmask = (dop>dopTh)&(bsxfun(@minus,(1:size(dop,1))',cath(3,:))>0);
%mmask(1:refInds(1)-1,:) = 0;

mmask = imdilate(imerode(mmask,ones(fwz,1)),ones(fwz,1));

dim = size(MM);
w = zeros([3,dim(2),dim(3)]);
if outputLevel>0
    w2 = zeros([3,dim(2),dim(3)]);
end

N = repmat([1;0;0;0;1;0;0;0;1],[1,1,dim(3)]);
%N = permute(makeRot(cat(1,zeros(2,numel(Vcorr)),-Vcorr)),[1,3,2]);

W = decomposeRot(MM);
Msqinv = makeRot(-W/2);
for indz = 2:roiz(end)%refInds(3):roiz(end)
    
    % D*N(indz)*D*Mtot(indz+1)*N'(indz)
    nloc = MatrixMultiply(bsxfun(@times,N,[1;1;-1;1;1;-1;-1;-1;1]),MatrixMultiply(MM(:,indz,:,:),N([1,4,7,2,5,8,3,6,9],:,:,:)));

    Wloc = decomposeRot(nloc);
 
    w(:,indz,:) = squeeze(Wloc);    
    nlocroot = makeRot(Wloc/2);

    Nnext = MatrixMultiply(nlocroot,N);
    N(:,:,mmask(indz,:)>0) = Nnext(:,:,mmask(indz,:)>0);

    if outputLevel>0
        % without depth correction
        nloc = MatrixMultiply(Msqinv(:,indz-1,:,:),MatrixMultiply(MM(:,indz,:,:),Msqinv(:,indz-1,:,:)));
%        nloc = MatrixMultiply(MM(:,indz,:,:),MM([1,4,7,2,5,8,3,6,9],indz-1,:,:));%
%        %this results in a optic axis with a V-component!

        Wloc = decomposeRot(nloc);
        w2(:,indz,:) = squeeze(Wloc);
    end
end

if outputLevel>0
    % without depth correction
    Omegaf = imfilter(permute(w2,[2,3,1]),ones(fwz,1)/fwz);
    phi2 = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));
    out.phi2 = phi2;
end
    
%Omega = bsxfun(@times,permute(pa,[2,3,1]),ret);
if ~cumulative % depth-resolved signal
    Omegaf = imfilter(permute(w,[2,3,1]),ones(fwz,1)/fwz);
    %retFinal = sqrt(sum(Omegaf.^2,3))/dzres/pi*180*100;
    phi = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));

    Omegafn = bsxfun(@rdivide,Omegaf,sqrt(sum(Omegaf.^2,3)));
else % outputting the cumulative signal
    Omegaf = permute(W,[2,3,1]);
    %retFinal = sqrt(sum(Omegaf.^2,3))/pi*100;% scale pi to 100
    phi2 = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));
    out.phi2 = mod(phi2 + pi,2*pi)-pi;
    Omegafn = bsxfun(@rdivide,Omegaf,sqrt(sum(Omegaf.^2,3)));
end
Omegaf = imfilter(permute(w,[2,3,1]),ones(fwz,1)/fwz);
%retFinal = sqrt(sum(Omegaf.^2,3))/dzres/pi*180*100;
phi = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));

out.dop = dop;
out.mask = mmask;
out.ret = retFinal;
out.phi = mod(phi + pi,2*pi)-pi;%do not include sheathAngle, as it may vary in between B-scans
%out.phi = mod(phi + out.sheathAngle + pi,2*pi)-pi;
out.PA = Omegafn;


out.fwz = fwz;
out.fwx = fwx;
out.unwrapMask = unwrapMask;
out.tom = tom;
out.MM = MM;

