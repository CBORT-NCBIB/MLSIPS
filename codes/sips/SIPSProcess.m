function out = SIPSProcess(S1,procStruct,Q,C)
% This script extimates the polarization properties of samples using a
% single input state
% INPUTS: Spectrally-binned Stokes vectors, system compensation matrices,
% processing parameters
% OUTPUTS: retardance, DOP, singular values to compute error metric,
% non-depth resolved optic axis

% the only two mandatory arguments
fwx = procStruct.fwx;
dz = procStruct.dz;
fwaxial = 1;% for comparison without spectral binning
fwz = 5; % axial filtering range of local ret
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

Omegaf = imfilter(permute(locoa,[2,3,1]),ones(dz,1)/dz);
retFinal = sqrt(sum(Omegaf.^2,3))/dzres/pi*180*100;

out = struct;
out.dop = sum(dop,3);
out.locoa = locoa;
out.ret = retFinal;
out.singVals = singVals;

end