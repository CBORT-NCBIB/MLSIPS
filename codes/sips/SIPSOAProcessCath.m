function out = SIPSOAProcessCath(S1,int,procStruct,Q,C)
% out = SIPSOAProcessCath(S1,int,procStruct,Q,C)
% Estimates the polarization properties of samples using a
% single input state, with depth-corrected optic axis and catheter
% transmission compensation

% INPUTS: Spectrally-binned Stokes vectors, intensity mask for catheter 
% system compensation matrices, processing parameters
% OUTPUTS: retardance, DOP, singular values to compute error metric,
% depth-resolved optic axis

% the only two mandatory arguments
fwx = procStruct.fwx;
fwz = 1; % default axial filtering range of local ret
dopTh = 0.7;
dzres = 4.8;% default axial sampling (um/px)

roiz = 1:1000; % axial range to consider for reconstruction of local ret
fw = 30; % lateral filtering to estimate surface signal and rotation due to sheath birefringence
cumulative = true; % flag for processing cumulative signal
outputLevel = 2;


fnames = fieldnames(procStruct);
for ind = 1:numel(fnames)
    if strcmp(fnames{ind},'fwz')
        fwz = procStruct.fwz;
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
h = h/sum(h(:));
S1f = imfilter(S1,h,'circular');

% final normalization
QUVf = sqrt(dot(S1f(:,:,:,2:4),S1f(:,:,:,2:4),4));
dop = mean(QUVf./S1f(:,:,:,1),3);

% Make Stokes into Jones vector for GDA
tom = Stokes2Jones(permute(S1f(:,:,:,2:4),[4,1,2,3]));

Nbins = size(S1,3);

dimQ = size(Q);

if dimQ(1) ~= 2
    C = reshape(C,[2,2,Nbins]);
end

r= svdDecomposition(C,Q,tom);
MM = squeeze(makeRot(r));

if outputLevel > 1
    out.OA1 = decomposeRot(MM);
end

% compute depth-resolved birefringence without optic axis orientation
dmn = MatrixMultiply(MM(:,2:end,:),MM([1,4,7,2,5,8,3,6,9],1:end-1,:));
locoa = cat(2,zeros(3,1,size(MM,3)),decomposeRot(dmn));
Omegaf = imfilter(permute(locoa,[2,3,1]),ones(fwz,1)/fwz);
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

% extract surface signal
mask = dop>0.8;
mask(1,:) = false; % make sure, at least the very first line is zero
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

% generate mask to exclude points with low DOP from correcting the local
% retardation
mmask = (dop>dopTh)&(bsxfun(@minus,(1:size(dop,1))',cath(3,:))>0);

mmask = imdilate(imerode(mmask,ones(fwz,1)),ones(fwz,1));

dim = size(MM);
w = zeros([3,dim(2),dim(3)]);
if outputLevel>0
    w2 = zeros([3,dim(2),dim(3)]);
end

N = repmat([1;0;0;0;1;0;0;0;1],[1,1,dim(3)]);

W = decomposeRot(MM);
Msqinv = makeRot(-W/2);
% iterative reconstruction with recovery of depth-resolved optic axis
% orientation
for indz = 2:roiz(end)
    
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
        Wloc = decomposeRot(nloc);
        w2(:,indz,:) = squeeze(Wloc);
    end
end

if outputLevel>0
    if ~cumulative
        % without depth correction
        Omegaf = imfilter(permute(w2,[2,3,1]),ones(fwz,1)/fwz);
        phi2 = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));
        out.phi2 = phi2;
    else % outputting the cumulative signal
        Omegaf = permute(W,[2,3,1]);
        phi2 = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));
        out.phi2 = mod(phi2 + pi,2*pi)-pi;
    end
end

Omegaf = imfilter(permute(w,[2,3,1]),ones(fwz,1)/fwz);
phi = atan2(real(Omegaf(:,:,2)),real(Omegaf(:,:,1)));

out.dop = dop;
out.mask = mmask;
out.ret = retFinal;
out.phi = mod(phi + pi,2*pi)-pi;




