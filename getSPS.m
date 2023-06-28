%% Characterize System Compensation
function[sps] = getSPS(sysComAv)
% Input: systemCompensation Structure
% Output: spectral polarization spread amount
% This function computes the overall spectral polarization spread present
% in the system. This can be used to see if the system has sufficiently 
% diverse polarization states interacting with the sample to accurately
% compute the poalrization properties


% Get rotation vector for illumination path (symmetric) system compensation
r = double(sysComAv.alignRotVec);

% Make rotation matrix
Mq = makeRot3x3(r);

% Find states at sample assuming horizontal input state
bp = squeeze(Mq(1,:,:));
% Find difference between adjacent bins
diffbp = bp(:,2:end)-bp(:,1:end-1);
% Get overall angular length
sps = sum(sqrt(sum(diffbp.^2,1)));
end