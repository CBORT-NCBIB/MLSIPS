function [r,singVals] = svdDecomposition(C,Q,eqqn)

Nbins = size(C,3);

C = reshape(C,[2,2,1,1,Nbins]);
Q = reshape(Q,[2,2,1,1,Nbins]);
eqqn = permute(eqqn,[1,5,2,3,4]);

t = pagemtimes(C,eqqn);
M = pagemtimes(conj(Q),[t(1,:,:,:,:).*conj(Q(1,1,:,:,:)),t(1,:,:,:,:).*conj(Q(1,2,:,:,:));...
                       t(2,:,:,:,:).*conj(Q(1,1,:,:,:)),t(2,:,:,:,:).*conj(Q(1,2,:,:,:))]);

dim = size(M);
M = reshape(M,[4,1,dim(3:end)]);

MM = -sum([M.*conj(M(1,:,:,:,:)),M.*conj(M(2,:,:,:,:)),M.*conj(M(3,:,:,:,:)),M.*conj(M(4,:,:,:,:))],5);

P = [1,-1i,0;0,0,-1i;0,0,-1i;1,1i,0]/sqrt(2);
H = pagemtimes(P',pagemtimes(MM,P));

% Eigenvalue decomposition
[U,S,V] = pagesvd(real(H));

Jopt = pagemtimes(P,-U(:,1,:,:))*sqrt(2);
r = JonesDecomp(Jopt);
singVals =cat(3,squeeze(S(1,1,:,:)),squeeze(S(2,2,:,:)),squeeze(S(3,3,:,:)));
end