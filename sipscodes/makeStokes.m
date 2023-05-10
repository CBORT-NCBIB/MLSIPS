function varargout = makeStokes(tom,outdim)
% t must be start with dim 2 
dimtom = [size(tom,2),size(tom,3),size(tom,4)];
dimtom(2) = dimtom(2)/2;

S0 = squeeze(abs(tom(1,:,:,:)).^2 + abs(tom(2,:,:,:)).^2);
S1 = squeeze(abs(tom(1,:,:,:)).^2 - abs(tom(2,:,:,:)).^2);
S2 = squeeze(2*real(tom(1,:,:,:).*conj(tom(2,:,:,:))));
S3 = squeeze(-2*imag(tom(1,:,:,:).*conj(tom(2,:,:,:))));

if nargout == 2
    
    varargout{1} = squeeze(reshape(cat(4,S0(:,1:2:2*dimtom(2),:),S1(:,1:2:2*dimtom(2),:),S2(:,1:2:2*dimtom(2),:),S3(:,1:2:2*dimtom(2),:)),[dimtom,4]));
    varargout{2} = squeeze(reshape(cat(4,S0(:,2:2:2*dimtom(2),:),S1(:,2:2:2*dimtom(2),:),S2(:,2:2:2*dimtom(2),:),S3(:,2:2:2*dimtom(2),:)),[dimtom,4]));

else
    if outdim == 3
    S = cat(4,S1,S2,S3);
    else
    S = cat(4,S0,S1,S2,S3);
    end
    varargout{1} = S;
end
