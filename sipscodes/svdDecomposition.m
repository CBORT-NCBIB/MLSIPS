function [r,singVals] = svdDecomposition(C,Q,eqqn)

Nbins = size(C,3);

for indw = 1:Nbins
    t(:,:,:,indw) = [C(1,1,indw).*eqqn(1,:,:,indw)+ C(1,2,indw).*eqqn(2,:,:,indw);C(2,1,indw).*eqqn(1,:,:,indw)+ C(2,2,indw).*eqqn(2,:,:,indw)];
end

for indw = 1:Nbins
    M(:,:,:,indw) = MatrixMultiply(conj(Q(:,indw)),[t(1,:,:,indw).*conj(Q(1,indw));t(2,:,:,indw).*conj(Q(1,indw));t(1,:,:,indw).*conj(Q(2,indw));t(2,:,:,indw).*conj(Q(2,indw))]);
end

conjM = conj(M);

for ind1 = 1:4
    for ind2 = 1:4
        MM(ind1,ind2,:,:) = -sum(permute(squeeze(M(ind1,:,:,:)),[3,1,2]).*permute(squeeze(conjM(ind2,:,:,:)),[3,1,2]));
    end
end

P = [1,-1i,0;0,0,-1i;0,0,-1i;1,1i,0]/sqrt(2);
H = pagemtimes(P',pagemtimes(MM,P));

% Eigenvalue decomposition
[U,S,V] = pagesvd(real(H));

Jopt = pagemtimes(P,-U(:,1,:,:))*sqrt(2);
r = JonesDecomp(Jopt);
singVals =cat(3,squeeze(S(1,1,:,:)),squeeze(S(2,2,:,:)),squeeze(S(3,3,:,:)));

end