function out = SIPSProcess(S1,procStruct,Q,C)
% out = SIPSOAProcessCath(S1,int,procStruct,Q,C)
% Estimates the polarization properties of samples using a single input 
% state
%
% INPUTS: Spectrally-binned Stokes vectors, system compensation matrices,
% processing parameters
% OUTPUTS: retardance, DOP, singular values to compute error metric,
% non-depth resolved optic axis

% the only two mandatory arguments
fwx = procStruct.fwx;
fwz = 1; % default axial filtering range of local ret vector
dzres = 4.8;

fnames = fieldnames(procStruct);
for ind = 1:numel(fnames)
    if strcmp(fnames{ind},'fwz')
        fwz = procStruct.fwz;
    elseif strcmp(fnames{ind},'dzres')
        dzres = procStruct.dzres;
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

%tom = tom./sqrt(abs(tom(1,:,:,:)).^2+abs(tom(2,:,:,:)).^2);
Nbins = size(S1,3);

dimQ = size(Q);

if dimQ(1) ~= 2
    C = reshape(C,[2,2,Nbins]);
end

[r,singVals]= svdDecomposition(C,Q,tom);

MM = squeeze(makeRot3x3(r));
MMT = pagetranspose(MM);
dmn = pagemtimes(MM(:,:,2:end,:),MMT(:,:,1:end-1,:));
locoa = cat(2,zeros(3,1,size(MM,4)),decomposeRot(dmn));

Omegaf = imfilter(permute(locoa,[2,3,1]),ones(fwz,1)/fwz);
retFinal = sqrt(sum(Omegaf.^2,3))/dzres/pi*180*100;%retFinal is in units of degrees of retardance per 100µm of tissue depth

out = struct;
out.dop = sum(dop,3);
out.locoa = locoa;
out.ret = retFinal;
out.singVals = singVals;

end