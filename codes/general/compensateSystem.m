function [MM,sysCom,MMcorr] = compensateSystem(MM,mask,sysComIn,computeErrBoolean)
%[MMcorr] = compensateSystem(MM,mask,sysCom) returns the symmetrized and 
%system-compensated SO3 rotation matrices, averaged across the spectral 
%bins.
%Size of MM is [3,3,Nz,NAlines,Nbins] or [9,Nz,NAlines,Nbins] and each 3x3 
%matrix is a pure SO3 rotation matrix. mask has a size of [Nz,NAlines] and
%defines the points to be used to compute the correction coefficients.
%The optional sysComIn is a systemCompensation object, defining the
%correction terms. computeErrBoolean triggers the computation of the errors
%even if correction terms are provided, to confirm if a given correction is
%still suitable.
%
%MM is the corrected output SO3 matrix, averaged across spectral bins.
%sysCom is the corresponding (updated) ystemCompensation object. MMcorr is
%the compenstated, but unaveraged SO3 matrices.
%
%
% Last modified October 2021 Martin Villiger


if nargin<3 || isempty(sysComIn)
    sysCom = systemCompensation;
else
    sysCom = sysComIn;
end

if nargin<4
    computeErrBoolean = false;
end

if ~isempty(MM)
    dim = size(MM);
    dimIn = dim;
    if dim(1) == 3 % backward compatibility; Here we want to use the first dimension to be 9
        MM = reshape(MM,[9,dim(3:end)]);
        dim = size(MM);
    end

    if numel(dim)<=3 % manage exception of no spectral binning
        dim(4) = 1;
    end
    
    if computeErrBoolean% force computation of error, even previous Hmat and alignProdSum exist
        sysCom.Hmat = [];
        sysCom.alignSumProd = []; 
    end
end

% Compute the correction necessary for making the system symmetric. Because
% this is quite fast, compute it anyway, even if wcorrIn is provided, for
% book-keeping over an entire volume. However, if wcorr and rcIn are
% provided, don't apply the correction, yet.

% Symmetrize matrices
if isempty(sysCom.symRotVec)% either sysCom is entirely empty or has an H-matrix
    [MM,sysCom] = makeSymmetric(MM,mask,sysCom);
elseif computeErrBoolean
    % don't apply correction, only compute errors
    [~,sysCom] = makeSymmetric(MM,mask,sysCom,true);
end



% The bin-alignment is more computationally costly. Only perform it if rcIn
% is not available.
if isempty(sysCom.alignRotVec)
%    if outputLevel < 1
    if nargout>2
        [MM,sysCom,MMcorr] = alignToCentralBin(MM,mask,sysCom);
    else
        [MM,sysCom] = alignToCentralBin(MM,mask,sysCom);
    end
else
    if computeErrBoolean
        [~,sysCom] = alignToCentralBin([],mask,sysCom);
    end
    if ~isempty(MM)
    % in this case, MM has not yet been symmetrized, and we now have to
    % apply to combined symmetry and alignment correction, and suppress
    % the v-components
    MMcorr = MM;
    for indw = 1:dim(4)
        Rsym = makeRot3x3(sysCom.symRotVec(:,indw));
        Ralign = makeRot3x3(sysCom.alignRotVec(:,indw));
        RalignT = Ralign.'.*[1,1,-1;1,1,-1;-1,-1,1];
        % What we need is RalignT*Rsym*MM*Ralign; we use A*B*C =
        % kron(C.',A)*vec(B)
        RR = kron(Ralign.',RalignT*Rsym);
        MMcorr(:,:,:,indw) = reshape(RR*reshape(MM(:,:,:,indw),[9,dim(2)*dim(3)]),[9,dim(2),dim(3)]);
    end
    
    MM = mean(MMcorr,4);
    MM = euclideanRotation(MM + MM([1,4,7,2,5,8,3,6,9]',:,:).*[1;1;-1;1;1;-1;-1;-1;1],true);

    %   Potentially a little bit faster, but less accurate...    
    %     mm = mean(decomposeRot(MMcorr),4);    
    %     oldphi = sqrt(sum(mm.^2,1));
    %     newphi = acos(cos(oldphi/2)./sqrt(1-sinc(oldphi/2/pi).^2.*mm(3,:,:,:).^2/4))*2;
    %     mm(3,:) = 0;
    %     mm = mm./sqrt(sum(mm.^2,1)).*newphi;
    %     MM2 = makeRot3x3(mm);

    %    outStruct.alignErrInit = [];
    %    outStruct.alignErrEff = squeeze(sum(sum(sum((MMcorr-reshape(MM,dim(1:3))).^2.*shiftdim(mask,-1),3),2),1))/sum(sum(mask));    
    %    outStruct.alignRotVec = [];
    %    outStruct.alignRotVecEff = rcIn;

    end
end

if ~isempty(MM)
    MM = reshape(MM,dimIn(1:end-1));
end