function [Mout,sysCom,Mcorr] = alignToCentralBin(Min,mask,sysComIn)
% [Mout,sysCom,Mcorr] = alignToCentralBin(Min,mask,sysComIn)
% alignToCentralBin estimates the required correction matrix to align all 
% spectral bins to the central bin. Only points that are true within mask 
% are used. Min is a 3 x 3 x Nz x NAlines x Nbins SO(3) symmetric (linear) 
% rotation matrix. The former 9 x Nz x NAlines x Nbins is also acceptable. 
%
% The optoinl input argument sysComIn is either a alignRotVecIn(3 x Nbins),
% defining the correction vector, to impose a defined alignment, or a 
% systemCompensation object. Optional output arguments are the (updated) 
% systemCompensation object, and Mcorr, which is the aligned SO3 matrices
% but without averaring among the spectral bins.
% If called with sysComIn containing a precomputed alignSumProd matrix, MM
% and mask can be empty, and a suitable alignRotVec is computed based on
% alignSumProd. If a alignRotVec is provided, but no alignSumProd, this
% will trigger computing a new alignProdSum and the corresponding error
% metrics, but use the provided alignRotVec.
%
%
% Detailed description:
% We are looking for an SO3 matrix L that minimizes the Frobenius norm of 
% D*L.'*D*Mbin*L-McentralBin, summed over all available masked points.
% Using the Kronecker product rule, this is equivalent to minimizing
% -trace(vec(kron(L.',D*L.'*D)).'*sum(vec(McentralBin)*vec(Mbin).')) which 
% is optimized here using fminsearch.
%
%
% Last modified October 2021 Martin Villiger


if nargin<3 || (~isobject(sysComIn) && isempty(alignRotVecIn))
    sysCom = systemCompensation;
elseif ~isobject(sysComIn)
    sysCom = systemCompensation;
    sysCom.alignRotVec = sysComIn;
else
    sysCom = sysComIn;
end


if isempty(sysCom.alignSumProd)
    % throughout this function, the first dimension will be forced to 9
    dim = size(Min);
    if dim(1) == 3 % for backward compatibility; convert 3x3 to linearized matrices
        Min = reshape(Min,[9,dim(3:end)]);
        dim = size(Min);
    end

    if numel(dim)<=3% manage exception of no spectral binning
        Mout = Min;
        Mcorr = [];
        return;
    end
    Nbins = dim(4);
    indwc = ceil(Nbins/2);     
    
    % apply mask to central bin
    McentralBin = reshape(Min(:,:,:,indwc),[9,dim(2)*dim(3)]);
    McentralBin = McentralBin(:,mask(:));
    Npoints = sum(mask(:));

    % compute alignSumProd
    alignSumProd = repmat(eye(9),[1,1,dim(4)]);
    for indw = [indwc+1:dim(4),indwc-1:-1:1]
        % apply mask to current bin
        Mbin = reshape(Min(:,:,:,indw),9,dim(2)*dim(3));
        Mbin = Mbin(:,mask);

        alignSumProd(:,:,indw) = (McentralBin*Mbin')/Npoints;
    end
    sysCom.alignSumProd = alignSumProd;
    sysCom.NSumProd = Npoints;
else
    alignSumProd = sysCom.alignSumProd;
    Nbins = size(alignSumProd,3);
    indwc = ceil(Nbins/2);     
end

Mcorr = Min;

LL = @(x)kron(makeRot3x3(x).',makeRot3x3(x).'.*[1,1,-1;1,1,-1;-1,-1,1]);
RR = repmat(eye(9),[1,1,Nbins]);
if isempty(sysCom.alignRotVec)
    xopt = [0;0;0];
    alignRotVec = zeros(3,Nbins);
    for indw = [indwc+1:Nbins,indwc-1:-1:1]
        fun = @(x)-trace(LL(x).'*alignSumProd(:,:,indw));

        xopt = fminsearch(fun,xopt);
        alignRotVec(:,indw) = xopt;
        if indw == Nbins
            xopt = [0;0;0];
        end

        % more efficient implementation of pre- and post-multiplying with
        % the same matrix
        RR(:,:,indw) = LL(alignRotVec(:,indw));
%        Mcorr(:,:,:,indw) = pagemtimes(RR,Min(:,:,:,indw));
    end
    sysCom.alignRotVec = alignRotVec;
else
    alignRotVec = sysCom.alignRotVec;
    for indw  = [indwc+1:Nbins,indwc-1:-1:1]
        RR(:,:,indw) = LL(alignRotVec(:,indw));
    end
end

if ~isempty(Min)
    for indw  = [indwc+1:Nbins,indwc-1:-1:1]
        Mcorr(:,:,:,indw) = pagemtimes(RR(:,:,indw),Min(:,:,:,indw));
    end
    % compute Euclidean mean rotation
    Mout = euclideanRotation(mean(Mcorr,4));
else
    Mout = [];
end

if nargout>1
    for indw = 1:Nbins
        % compute the errors, which correspond to 6 -
        % 2*trace(LL.'*alginSumProd)
        sysCom.alignErrInit(indw) = 6-2*trace(alignSumProd(:,:,indw));
        fun = @(x)-trace(LL(x).'*alignSumProd(:,:,indw));
        sysCom.alignErr(indw) = 6+2*fun(sysCom.alignRotVec(:,indw));
    end
    sysCom.alignErrInit(indwc) = 0;
    sysCom.alignErr(indwc) = 0;
end

