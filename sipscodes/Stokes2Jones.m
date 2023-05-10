function jones = Stokes2Jones(stokes)
%tom = Stokes2Jones(stokes) transforms a Stokes vector into a matching Jones
%vector. The input Stokes vectors should be a three components Stokes
%vector aligned along the first dimension.

% The conversion follows the following convention, where h and v are the
% two components of the Jones vector.
% 
% Q = |h|^2 - |v|^2
% U = 2*real(h*conj(v))
% V = -2*imag(h*conj(v))
%
% Hence, I = sqrt(Q^2 + U^2 + V^2) = |h|^2 + |v|^2, and
% |h| = sqrt((I + Q)/2), and
% U + 1i*V = 2*|h|*|v|*exp(1i*phi), with phi angle(v/h)
% from this: |v| = (U+1i*V)/2/|h|

dim = size(stokes);
if dim(1)~=3
    error('Input data must have three components along the first dimension');
end

I = sqrt(sum(stokes(:,:).^2,1));
a = sqrt((I+stokes(1,:))./2);
b = (stokes(2,:) + 1i*stokes(3,:))/2./a;

% reshape the output Jones vectors according to the other dimensions of the
% input Stokes vector
jones = reshape(cat(1,a,b),[2,dim(2:end)]);
    