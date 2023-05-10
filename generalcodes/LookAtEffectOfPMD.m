%% Look at effect of PMD

%% File & Parameters
fname = '/Volumes/LaCie SSD/ClinicalData/142988_12-16-2014_04-55-47';
startSlice = 100;
numSlices = 10;
sliceInds = startSlice:startSlice+(numSlices-1);
st.window = 5;
Nbins = 2*st.window-1;
st.skipLastPoints = 30; % this is an adhoc fix that helps to suppress artifacts in some measurements
st.skip = 1;
noiseLevel = 0.04;


%% Get QC Estimation
numSlicesComputation = 20;
%sysCompSIPS = estimateQCSIPS(fname,numSlicesComputation,startSlice,Nbins)

QSIPS = reshape(makeJones(-sysCompSIPS.alignRotVec),[2,2,Nbins])
CSIPS = reshape(makeJones(sysCompSIPS.symRotVec),[2,2,Nbins])

%%
foutPMD = "WithPMD.pstif";
foutNoPMD = "NoPMD.pstif";
maxRet = 200;
logLim = [55,110];
fwx = 12;
fwz = [];
dz = 5;
N = 9;

fid = fopen(foutPMD,'w+');
fclose(fid);
%%
for ind = 1:size(sliceInds,2)
disp("DIPS")
[S1w,S2w] = recstrTom(fname,sliceInds(ind),st);
[S1,S2]= recstrTom(fname,sliceInds(ind),struct);

int = 10*log10(tom2Int(S1,S2));

dz = 5;
dzres = 4.8;
pstruct = struct;
pstruct.fwx = 12; 
pstruct.dz = 5;
pstruct.dzres = 4.8;

[out,OA] = PSProcessGlobalSymmetric(S1w,S2w,pstruct);

%% Get cumulative retardance matrix
MM = out.MM;

% To Jones
MMJ = makeJones(decomposeRot(MM));

% Add Noise
MMJNoise = MMJ + noiseLevel*(rand(size(MMJ))+1j*rand(size(MMJ)));

% Get First Column
tCorr = MMJNoise(1:2,:,:);

% Complete Jones Matrix
JQQ = permute(cat(3,squeeze(tCorr(1,:,:)),squeeze(tCorr(2,:,:)),squeeze(tCorr(2,:,:)),squeeze(-conj(tCorr(1,:,:))).*squeeze(exp(2*1i*angle(tCorr(2,:,:))))),[3,1,2]);

MMCorr = reshape(makeRot(JonesDecomp(JQQ)),[3,3,size(JQQ,2),size(JQQ,3)]);

dmn = pagemtimes(MMCorr(:,:,2:end,:),pagetranspose(MMCorr(:,:,1:end-1,:)));
locoa = cat(2,zeros(3,1,size(MM,4)),decomposeRot(dmn));

Omegaf = imfilter(permute(locoa,[2,3,1]),ones(dz,1)/dz);
retFinal = sqrt(sum(Omegaf.^2,3))/dzres/pi*180*100;

%% Plot
dopThresh = 0.7;
retOG = out.ret;
retNoPMD = retFinal;

retOG(out.dop<dopThresh) = 0;
retNoPMD(out.dop<dopThresh) = 0;

figure,
subplot(1,2,1)
imagesc(retOG,[0,120])
axis image
colorbar
title("With PMD")
subplot(1,2,2)
imagesc(retNoPMD,[0,120])
axis image
colorbar
title("No PMD")

%% Write as pstif

% 
% format = 'IRD';
% Format = {'Int:',logLim(1),logLim(2),[];'Ret:',0,maxRet,[];'DOP:',0,1,[];'ParamsNFwxDz',N,fwx,dz};
% tag = sprintf('Format:\t%s\n',format);
% for ind = 1:numel(format)
%     tag = [tag,sprintf('%s\t%d\t%d\n',Format{ind,1},Format{ind,2},Format{ind,3})];
% end
% intout = uint8(255*(int-logLim(1))/diff(logLim));
% 
% %PMD
% retout = uint8(255*out.ret/maxRet);
% dopout = uint8(255*out.dop);
% imwrite(intout,foutPMD,'tif','WriteMode','append','Compression','none','Description',tag);
% imwrite(retout,foutPMD,'tif','WriteMode','append','Compression','none','Description',tag);
% imwrite(dopout,foutPMD,'tif','WriteMode','append','Compression','none','Description',tag);
% 
% % No PMD
% retout = uint8(255*retFinal/maxRet);
% dopout = uint8(255*out.dop);
% imwrite(intout,foutNoPMD,'tif','WriteMode','append','Compression','none','Description',tag);
% imwrite(retout,foutNoPMD,'tif','WriteMode','append','Compression','none','Description',tag);
% imwrite(dopout,foutNoPMD,'tif','WriteMode','append','Compression','none','Description',tag);
end