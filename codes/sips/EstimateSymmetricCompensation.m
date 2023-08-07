

function qOut =  EstimateSymmetricCompensation(CT,N,S1f)

dopThresh = 0.8;

%[St,dop] = getS1(fname,startSlice);
[St,dop] = getDOP(S1f);
St = St(:,:,:,2:4);

% Apply C to S1 Vector
for i = 1:N
    CTbin = squeeze(CT(:,:,i));
    Stbin = squeeze(St(:,:,i,:));
    SCorr(:,:,i,:) = permute(pagemtimes(CTbin,permute(Stbin,[3,1,2])),[2,3,1]);
end

% Convert to Jones 
tCorr = Stokes2Jones(permute(SCorr,[4,1,2,3]));

% Get JQQ
JQQ = cat(4,squeeze(tCorr(1,:,:,:)),squeeze(tCorr(2,:,:,:)),squeeze(tCorr(2,:,:,:)),squeeze(-conj(tCorr(1,:,:,:))).*squeeze(exp(2*1i*angle(tCorr(2,:,:,:)))));

% Compute MQQ (back to stokes)
rQQ = JonesDecomp(permute(JQQ,[4,1,2,3]));
MQQ = makeRot(rQQ);

[Mout,sysCompensation] = alignToCentralBin(MQQ,dop>dopThresh);%findAlignMatrix(MQQ,dop,0.6);

% Initialize Alignment Matrix & Vector
QSIPS = makeRot(sysCompensation.alignRotVec);
QvecSIPS = sysCompensation.alignRotVec;

qOut.Q = QSIPS;
qOut.qVec = QvecSIPS;
qOut.sysCompensation = sysCompensation;
qOut.Mout = Mout;
qOut.MQQ = MQQ;
qOut.mask = dop>dopThresh;

end
function [S1n,dop] = getS1(fname,sliceInds)   
fwx = 12;
dz = 5;
fwz = 1;
dopTh = 0.6;
wcorr = [];
rc = [];
dzres = 4.8;

for i = 1:size(sliceInds,2)
    % Get S1
    sliceInd = sliceInds(i);
    %logF = readLogFile(fname);
    st = struct;
    st.skipLastPoints = 30;
    st.window = 5;
    [S1,S2] = recstrTom(fname,sliceInd,st);
    int = tom2Int(S1,S2);
    %% Filter & Get DOP
    dim = size(S1);

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
    %
    % final normalization
    QUVf = sqrt(dot(S1f(:,:,:,2:4),S1f(:,:,:,2:4),4));
    dop = mean(QUVf./S1f(:,:,:,1),3);

    S1n = (S1f./QUVf); 
end
end

function[S1n,dop] = getDOP(S1f)
    QUVf = sqrt(dot(S1f(:,:,:,2:4),S1f(:,:,:,2:4),4));
    dop = mean(QUVf./S1f(:,:,:,1),3);
    S1n = (S1f./QUVf); 
end
