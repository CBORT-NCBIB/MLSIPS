function [MMcorr,sysCom,MMcorrraw] = makeSymmetric(MM,mask,sysComIn,boolNoMMcorr)
%[MMcorr,sysCom,MMcorrraw] = makeSymmetric(MM,mask,wcorr) returns the 
%symmetrized measurement matrix. Size of MM is [9,Nz,NAlines,Nbins] or 
%[3,3,Nz,NAlines,Nbins] and each 3x3 matrix is a pure rotation matrix. mask 
%has a size of [Nz,NAlines]. The optional input argument sysComIn, which 
%can be either the rotationVector (for backward compatibility), or a 
%systemCompensation object results in directly applying the defined 
%correction, or recomputing the rotation vectors based on the provided Hmat 
%matrix.
% 
%MMcorr is the symmetrized MM, with the remaining component along the
%V-coordinate set to zero. If boolNoMMCorr is set to true, its computation
%is suppressed.
%If a sysComIn contains an Hmat field, MM and mask can be empty, simply
%prompting the computation of the corresponding symRotVec. If a symRotVec 
% is provided, but no Hmat, thiswill trigger computing a new Hmat and the 
% corresponding error metrics, but use the provided symRotVec.
%
%MMcorr is the symmetrized SO3 matrices, in 3x3xNzxNAlinesxNBins format. 
%sysCom is a systemCompensation object. MMcorrraw is the symmetrized MM, 
%but without setting the V-coordinate to zero.
%
%
% Last modified October 2021 Martin Villiger


if nargin<3 || ~isobject(sysComIn)
    sysCom = systemCompensation;
elseif ~isobject(sysComIn) % backward compatibility
    sysCom = systemCompensation;
    sysCom.symRotVec = sysComIn;
elseif isobject(sysComIn)
    sysCom = sysComIn;
else 
    return;
end

if nargin<4 || isempty(boolNoMMcorr)
    boolNoMMcorr = false;
end



if isempty(sysCom.Hmat)
    dim = size(MM);
    if dim(1) == 9 % backward compatibility
        MM = reshape(MM,[3,3,dim(2:end)]);
        dim = size(MM);
    end

    if numel(dim)<=4 % manage exception of no spectral binning
        dim(5) = 1;
    end

    % First, average the SO3 matrices over the mask ww
    Npoints = sum(mask(:));% number of points in mask
    MM(isnan(MM))=0;% make sure there are no nans present
    for jnd = 1:dim(1)
        for ind = 1:dim(2)
            %size(MMproj) is 9 by number of bins
            MMproj(ind+(jnd-1)*3,:) = mask(:)'*reshape(MM(ind,jnd,:,:,:),numel(mask),dim(5))/Npoints;
        end
    end

    % compute H matrix, unless already available
    A = [1 0 0 1; 1 0 0 -1; 0 1 1 0; 0 1i -1i 0]/sqrt(2);
    for indw = 1:dim(5)
        % find estimation matrix, convert MMproj to H-matrix, then multiply
        % with diag([-1,1,-1,1]) on both sides
        F = conj(A'*[1,0,0,0;cat(2,zeros(3,1),reshape(MMproj(:,indw),[3,3]))]*A);
%        F = conj(A'*[sum(mask(:)),0,0,0;cat(2,zeros(3,1),reshape(MMproj(:,indw),[3,3]))]*A);
        H(:,:,indw) = reshape(F([1,9,3,11,5,13,7,15,2,10,4,12,6,14,8,16]'),[4,4]); 
    end
    sysCom.Hmat = H;
    sysCom.NH = Npoints;
    Nbins = dim(5);
else
    H = sysCom.Hmat;
    Nbins = size(H,3);
end

% compute symRotVec
if isempty(sysCom.symRotVec)
    for indw = 1:Nbins
        [a,b] = eig(diag([-1,1,-1,1])*H(:,:,indw)*diag([-1,1,-1,1]));
        [~,pos] = min(abs(diag(b)));
        cc = a(:,pos);

        Jcorr = [cc(2),cc(4);cc(1),cc(3)];
        Jcorr = Jcorr/sqrt(det(Jcorr));

        rr = JonesDecomp(Jcorr);% convert back to Stokes, ignoring diattenuation component

        % collect correction rotation vector for all spectral bins
        symRotVec(:,indw) = real(rr);
%        ddcorr(:,indw) = real(dd);
    end

    % unwrap correction across spectral bins
    ss = cat(2,0,mod(cumsum(sqrt(sum(diff(symRotVec,[],2).^2,1))>pi),2));
    if sum(ss)>Nbins/2
        ss = 1-ss;
    end
    retcorr = sqrt(sum(symRotVec.^2,1));
    symRotVec = symRotVec - symRotVec./retcorr.*ss*2*pi;

    % maintain a retardation <pi
    if mean(sqrt(sum(symRotVec.^2,1)))>pi
        symRotVec = symRotVec - symRotVec./sqrt(sum(symRotVec.^2,1))*2*pi;
    end

    sysCom.symRotVec = symRotVec;
else
    symRotVec = sysCom.symRotVec;
end

% Apply correction; only do this, if MM is requested as output argument
if ~boolNoMMcorr && ~isempty(MM)
    % attempt to speed this up
    Rcorr = reshape(makeRot3x3(symRotVec),[3,3,1,1,Nbins]);
    MMcorr = pagemtimes(Rcorr,MM);

    if nargout>2
        % maintain MMcorr in memory before setting V-components to zero 
        MMcorrraw = MMcorr;
    end

    % cancel out V-component precisely, by scaling the overall retardance
    % to take into account the clipped V-component. This is equivalent to
    % taking the coherent mean of J + J.' in Jones, or euclideanRotation(M
    % M+D*M.'*D) in SO3
    OA = decomposeRot(MMcorr);
    oldphi = sqrt(sum(OA.^2,1));
    newphi = real(acos(cos(oldphi/2)./sqrt(1-sinc(oldphi/2/pi).^2.*OA(3,:,:,:).^2/4))*2);
    OA(3,:) = 0;
    OA = OA./sqrt(sum(OA.^2,1)).*newphi;
    MMcorr = makeRot3x3(OA);
    %MMcorr = euclideanRotation(MMcorr + MMcorr([1,4,7,2,5,8,3,6,9]',:,:,:).*[1;1;-1;1;1;-1;-1;-1;1],true);
    
else
    % if boolNoMMcorr is true, simply assign empty matrix
    MMcorr = [];
    MMcorrraw = [];
end

if nargout>1
   for indw = 1:Nbins
        Jcorr = makeJones(symRotVec(:,indw));
        jvec = [Jcorr(2,1),Jcorr(1,1),Jcorr(2,2),Jcorr(1,2)].';
        errOpt(indw) = real(jvec'*diag([-1,1,-1,1])*H(:,:,indw)*diag([-1,1,-1,1])*jvec/4); % factor of 4 is as reference for entire energy; each pixel within ww is a normalized Jones matrix with sum(abs(J(:)).^2) = 2
        % more pragmatically, a pure V-component matrix M with pi
        % retardation results in an error of 1
        jvec = [0, 1, 1, 0]';
        errInit(indw) = real(jvec'*diag([-1,1,-1,1])*H(:,:,indw)*diag([-1,1,-1,1])*jvec/4);
    end
    sysCom.symErrInit = errInit;
    sysCom.symErr = errOpt;
end

