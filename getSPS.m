%% Characterize System Compensation
function[sps] = getSPS(sysComAv)

% Get rotation vector for illumination path (symmetric) system compensation
r = double(sysComAv.alignRotVec);

% Find SPS
Mq = makeRot3x3(r);
bp = squeeze(Mq(1,:,:));
diffbp = bp(:,2:end)-bp(:,1:end-1);
sps = sum(sqrt(sum(diffbp.^2,1)));
end