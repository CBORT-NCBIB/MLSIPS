


function [HH, varargout] = poincareHistogram(S,Naz,Nel,mask)
%% Poincare shpere digitaliation
binsAz = linspace(-pi,pi,Naz);
binsEl = linspace(-pi/2,pi/2,Nel);

elAngle = real(asin(S(:,:,:,3)));
azPos = real(atan2(S(:,:,:,2),S(:,:,:,1)).*cos(elAngle));

% figure(32);clf;
% subplot(3,2,1);imagesc(S(:,:,5,1).*mask);colorbar;
% subplot(3,2,3);imagesc(S(:,:,5,2).*mask);colorbar;
% subplot(3,2,5);imagesc(S(:,:,5,3).*mask);colorbar;
% subplot(3,2,2);histogram(S(:,:,5,1).*mask);xlim([-1,1]);ylim([0,1500])
% subplot(3,2,4);histogram(S(:,:,5,2).*mask);xlim([-1,1]);ylim([0,1500])
% subplot(3,2,6);histogram(S(:,:,5,3).*mask);xlim([-1,1]);ylim([0,1500])
% 
% Q=S(:,:,6,1);
% QM=Q(mask);
% U=S(:,:,6,2);
% UM=U(mask);
% V=S(:,:,6,3);
% VM=V(mask);
% 
% poincareQY(QM(1:3000),UM(1:3000),VM(1:3000),0);

HH=zeros(1001,2001,size(S,3));

for wind = 1:size(S,3)
    elLoc = real(elAngle(:,:,wind));
    azLoc = real(azPos(:,:,wind));
    hh = hist2D(elLoc(mask),azLoc(mask),binsEl,binsAz);
    hh = imfilter(hh,ones(2,4)/8,'symmetric');
    HH(:,:,wind) = hh;
end

if nargout>1
    varargout{1} = binsAz;
end
if nargout>2
    varargout{2} = binsEl;
end

