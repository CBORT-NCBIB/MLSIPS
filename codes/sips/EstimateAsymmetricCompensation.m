function[out,HHfinal] = EstimateAsymmetricCompensation(S1,Naz,Nel)
%%


binsAz = linspace(-pi,pi,Naz);
binsEl = linspace(-pi/2,pi/2,Nel);
Npoints = 20;

% Get histogram (HFinal) and sum of all S1 vectors
disp("Computing Histograms...")
[HHfinal] = getHistogramData(S1);

% % Visualize Histogram
% for i = 1:size(HHfinal,3)
%     figure('Renderer', 'painters','Position', [200,200 400 600]),
%     hBin = HHfinal(:,:,i);
%     % Find max
%     [maxVal,indsMax] = max(hBin,[],'all');
%     hBinNorm = hBin./maxVal;
%     [indsx,indsy]= ind2sub(size(hBin),indsMax);
%     % Find histogram points of MS
%     imagesc(hBinNorm,[0,1]),colormap('hot')
%     hold on
%     plot(indsy,indsx,'b*','MarkerSize',20)
%     set(gca,'LineWidth',1.5,'FontSize',15)
%     set(gca,'xtick',[20,size(binsAz,2)-20],'XTickLabel',["-\pi/2","\pi/2"],'ytick',[20,size(binsEl,2)-20],'YTickLabel',["-\pi","\pi"])
%     ylabel('\theta cos(\phi)')
%     xlabel('\phi ')
%     colorbar
%     %print("PSOFDI_Histograms_200_Bin"+i,'-djpeg','-r500')
% end

% Implement Randsac
[symMatrix,minErr] = randsac(HHfinal,binsAz,binsEl,Npoints);

out.symMatrix = symMatrix;
out.minError = minErr;

end

function[RR,currMinErr] = randsac(HHFinal,binsAz,binsEl,Npoints)
%%
N = 9;

% Initial State for R
rdcInit = rand(3,1);
rlinInit= rand(3,1);

% Initial vector which will be optimized
vecInit = [rdcInit;rlinInit];

params(1) = 10; % array size
params(2) = 0.6; % maxDifference
%% Find maxima
[msMax,numMax] = findMaximum2DHistogram(HHFinal,N,params,Npoints);

% Define search paramaters
disp("Using LSQ to compute C...")
maxIter = 200;
currMinErr = 10;
errTot = [];
RRTot = [];
counter = 0;
options = optimoptions(@lsqnonlin,'Display','none','FunctionTolerance',1e-3);
options.Algorithm = 'levenberg-marquardt';

for ind = 1:maxIter
   
    % Choose random maxima
    for nBin = 1:N
        % Get 2D ind
        maxNum = randi([1 numMax(nBin)]);
        maxCell = msMax{nBin};
        indMax = maxCell(maxNum,1:2);
        % Compute to Stokes
        azAngle = binsAz(indMax(2));
        elAngle = binsEl(indMax(1));
        S3 = sin(elAngle);
        azEff = azAngle / cos(elAngle);
        S1 = cos(azEff)*cos(elAngle);
        S2 = sin(azEff)*cos(elAngle);
        tempMS(:,nBin) = [S1,S2,S3]';
    end
    tempMS;
    % Define LSQ error function
    func = @(vec) lsqFunc(vec,tempMS);

    % Find ideal vector
    [xFinal,err] = lsqnonlin(func,vecInit,[],[],options);

     
    if err < currMinErr
        disp("ind = "+num2str(ind)+", err = "+num2str(err))
        currMS = tempMS;
        currXFinal = xFinal;
        currMinErr = err;
        counter = counter+1;
        errTot = [errTot,currMinErr];
        rdcFinal = currXFinal(1:3);
        rlinFinal = currXFinal(4:6);
        RDC = makeRot(rdcFinal);%reshape(makeRot(rdc),[3,3]);
        Rlin = makeRot(rlinFinal*linspace(0,1,N));%reshape(makeRot(rlin*linspace(0,1,N)),[3,3,N]);
        RR = MatrixMultiply(RDC,Rlin);%pagemtimes(RDC,Rlin);
        RRTot = cat(3,RRTot,RR);
    end

end
%%
% Find "Mirror State" of estimated C matrix
rdcFinal = currXFinal(1:3);
rlinFinal = currXFinal(4:6);
RDC = makeRot3x3(rdcFinal);%reshape(makeRot(rdc),[3,3]);
Rlin = makeRot3x3(rlinFinal*linspace(0,1,N));%reshape(makeRot(rlin*linspace(0,1,N)),[3,3,N]);
RR = pagemtimes(RDC,Rlin);
pointsModel = RR(1,:,:);

end

%% Extra Functions
function[msMax,numMax] = findMaximum2DHistogram(HHFinal,N,params,Npoints)

arraySize = params(1);
minDiff = params(2);
dim = size(HHFinal);
msMax = [];
numMax = [];
%% Cycle through bins
for nBin = 1:N
    msMaxBin = [];
    hBin = squeeze(HHFinal(:,:,nBin));

    % Normalize
    hBin = hBin ./ max(max(hBin));
    hBinFilt = imgaussfilt(hBin,1);
    [hGradMag,hGradDir]= imgradient(imgaussfilt(hBin,2));

    
    % Find max Npoints of filtered hBin
    [maxGrad, indMax] = maxk(hBinFilt(:),Npoints);
    [indsx,indsy]= ind2sub(size(hBinFilt),indMax);

    % Remove close together regions
    %%
    barrier = 10; % must be at least 5 pixels apart
    i = 1;
    complete = 0;
    while complete == 0
        toRemove = zeros(1,size(indsx,1));
        for j = i+1:size(indsx,1);         
            xdiff = abs(indsx(i)-indsx(j));
            ydiff = abs(indsy(i)-indsy(j));
            if(xdiff<barrier && ydiff<barrier)
                toRemove(j) = 1;
            else
                toRemove(j) = 0;
            end
        end

        indsx(toRemove==1) = [];
        indsy(toRemove==1) = [];
        if((i+1)>size(indsx,1))
            complete = 1;
        else
            i = i+1;
        end
    end
% figure,
% imagesc(hBinFilt)
% hold on
% plot(indsy,indsx,'sqr','MarkerSize',16,'LineWidth',2,'Color','r')
% pause

msMaxBin = [indsx,indsy];
msMax{nBin} = msMaxBin;
numMax(nBin) = size(indsx,1);
end

end

function [HHSum] = getHistogramData(S1ftot)   
n_slices = size(S1ftot,5);
Nel = 1001;
Naz = 2001;
HHSum = zeros(Nel,Naz,9);
    
for sliceInd = 1:n_slices
    S1f = squeeze(S1ftot(:,:,:,:,sliceInd));
    % final normalization
    QUVf = sqrt(dot(S1f(:,:,:,2:4),S1f(:,:,:,2:4),4));
    S1n = (S1f./QUVf);
    dop = mean(QUVf./S1f(:,:,:,1),3);
    mask = dop>0.7;

    % Get Histogram
    [HH,binsEl,binsAz] = poincareHistogram(squeeze(S1n(:,:,:,2:4)),Naz,Nel,mask);
    HHSum = HHSum+HH;
end
end

function[err] = lsqFunc(vec,points)
N = 9;
rdc = vec(1:3);
rlin = vec(4:6);
 RDC = makeRot(rdc);%reshape(makeRot(rdc),[3,3]);
 Rlin = makeRot(rlin*linspace(0,1,N));%reshape(makeRot(rlin*linspace(0,1,N)),[3,3,N]);
 RR = MatrixMultiply(RDC,Rlin);%pagemtimes(RDC,Rlin);
 
 pointsModel = RR(1:3,:);

 % Compute LSQ
 err = 0;
 for i = 1:N
    binErr = sum(abs(points(:,i)-pointsModel(:,i)).^2);
    err = err + binErr;
 end

 
end

function curr = get2DHistInd(p)
    Nel = 1001;
    Naz = 2001;
    S = p;
    elAngle = real(asin(S(3)));
    azPos = real(atan2(S(2),S(1)).*cos(elAngle));

    binsAz = linspace(-pi,pi,Naz);
    binsEl = linspace(-pi/2,pi/2,Nel);

    % Find Index
    [~,ind1] = min(abs(binsEl-elAngle));
    [~,ind2] = min(abs(binsAz-azPos));

    curr = [ind1;ind2];
   
end

