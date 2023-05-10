
function[out,HHfinal] = EstimateAsymmetricComponent(S1)
%%

binsAz = linspace(-pi,pi,Naz);
binsEl = linspace(-pi/2,pi/2,Nel);
Npoints = 20;

% Get histogram (HFinal) and sum of all S1 vectors
disp("Computing Histograms...")
[HHfinal,S1fTot] = getHistogramData(fname,Nel,Naz,sliceInds,pstruct);

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
out.S1f = S1fTot;

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

red = [1, 0, 0];
blue = [0,0,1];
numVals = N;
colorGrad = [linspace(red(1),blue(1),numVals)', linspace(red(2),blue(2),numVals)', linspace(red(3),blue(3),numVals)'];
% Define search paramaters
disp("Using LSQ to compute C...")
maxIter = 200;
currMinErr = 10;
errTot = [];
RRTot = [];
counter = 0;
options = optimoptions(@lsqnonlin,'Display','none','FunctionTolerance',1e-3);
options.Algorithm = 'levenberg-marquardt';
% Display
% figure,
% hold on
% for binInd = 1:N
%     maxCell = msMax{binInd};
%     for maxInd = 1:numMax(binInd)
%         indMax = maxCell(maxInd,1:2);
%         azAngle = binsAz(indMax(2));
%         elAngle = binsEl(indMax(1));
%         S3 = sin(elAngle);
%         azEff = azAngle / cos(elAngle);
%         S1 = cos(azEff)*cos(elAngle);
%         S2 = sin(azEff)*cos(elAngle);
% 
%         plot3(S1,S2,S3,'.','Color',colorGrad(binInd,:),'MarkerSize',20)
%     end
% end
for i = 1:size(HHFinal,3)
    maxCell = msMax{i};
    figure('Renderer', 'painters','Position', [200,200 400 600]), 
    hBin = HHFinal(:,:,i);
    [maxVal,indsMax] = max(hBin,[],'all');
    hBinNorm = hBin./maxVal;
    % Find histogram points of MS
    imagesc(hBinNorm,[0,1]),colormap('hot')
    hold on
    for maxInd = 1:numMax(i)
     indMax = maxCell(maxInd,1:2);    
    plot(indMax(2),indMax(1),'bo','MarkerSize',20)
    end
     set(gca,'LineWidth',1.5,'FontSize',15)
    set(gca,'xtick',[20,size(binsAz,2)-20],'XTickLabel',["-\pi/2","\pi/2"],'ytick',[20,size(binsEl,2)-20],'YTickLabel',["-\pi","\pi"])
    ylabel('\theta cos(\phi)')
    xlabel('\phi ')
    colorbar
%     %print("PSOFDI_Histograms_200_Bin"+i,'-djpeg','-r500')
 end
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
%         pointsModel = RR(1:3,:);
%         if counter < 7
%             plot3(pointsModel(1,:),pointsModel(2,:),pointsModel(3,:),'k.-','MarkerSize',20,'LineWidth',2)
%             %plot3((1,:),pointsModel(2,:),pointsModel(3,:),'k.-','MarkerSize',20,'LineWidth',2)
%         end
        RRTot = cat(3,RRTot,RR);
    end

end
% plot3(pointsModel(1,:),pointsModel(2,:),pointsModel(3,:),'.-','Color',[0,0,0],'MarkerSize',20,'LineWidth',2)   
% % use surf function to plot
% [x,y,z] = sphere;
% mesh(x,y,z,'FaceColor','none','EdgeColor',[0.3,0.3,0.3])
% xlabel('Q','FontSize',20,'Color','k')
% ylabel('U','FontSize',20,'Color','k')
% zlabel('V','FontSize',20,'Color','k')
% set(gca,'xtick',[-1,1],'ytick',[-1,1],'ztick',[-1,1]);

%%
% Find "Mirror State" of estimated C matrix
rdcFinal = currXFinal(1:3);
rlinFinal = currXFinal(4:6);
RDC = makeRot(rdcFinal);%reshape(makeRot(rdc),[3,3]);
Rlin = makeRot(rlinFinal*linspace(0,1,N));%reshape(makeRot(rlin*linspace(0,1,N)),[3,3,N]);
RR = MatrixMultiply(RDC,Rlin);%pagemtimes(RDC,Rlin);
pointsModel = RR(1:3,:);

 %Plot Progress
% figure,
% plot3(currMS(1,:),currMS(2,:),currMS(3,:),'r.-','MarkerSize',20,'LineWidth',2)
% hold on
% plot3(pointsModel(1,:),pointsModel(2,:),pointsModel(3,:),'b.-','MarkerSize',20,'LineWidth',2)
% alpha = 0.4;
% [x,y,z]=sphere(64);
% use surf function to plot
% a = 1;
% hSurface=surf(a*x,a*y,a*z);
% set(hSurface,'FaceColor',[1 1 1], ...
%   'FaceAlpha',alpha,'FaceLighting','gouraud','EdgeColor','none')
% camlight
% set(gca,'CameraPosition',[10.8068   0.1271    4.9838])
% axis equal
% axis tight
% grid off
% axis off
% aa = 1.0;
% plot3([0,aa],[0,0],[0,0],'k','LineWidth',10)
% plot3([0,0],[0,aa],[0,0],'k','LineWidth',10)
% plot3([0,0],[0,0],[0,aa],'k','LineWidth',10)
% legend("MS Computed C","Estimated C using LSQ")
% %%
% for i = 1:size(HHFinal,3)
% figure('Renderer', 'painters','Position', [100 100 600 1000]),
% hBin = HHFinal(:,:,i);
% % Find max
% [maxVal,indsMax] = max(hBin,[],'all');
% hBinNorm = hBin./maxVal;
% [indsx,indsy]= ind2sub(size(hBin),indsMax);
% indsLSQ = get2DHistInd(currMS(:,i));
% % Find histogram points of MS
% imagesc(hBinNorm,[0,1])
% colormap('hot')
% hold on
% plot(indsy,indsx,'sqr','MarkerSize',16,'LineWidth',2,'Color','r')
% plot(indsLSQ(2),indsLSQ(1),'sqr','MarkerSize',16,'LineWidth',2,'Color','g')
% %colorbar;
% set(gca,'LineWidth',1.5,'FontSize',15)
% set(gca,'xtick',[],'ytick',[]);
% %legend("Maximum Points","LSQ Predicted")
% %title("Bin:"+num2str(i))
% set(gca,'LineWidth',1.5,'FontSize',15)
% set(gca,'xtick',[20,size(binsAz,2)-20],'XTickLabel',["-\pi/2","\pi/2"],'ytick',[20,size(binsEl,2)-20],'YTickLabel',["-\pi","\pi"])
% ylabel('\theta cos(\phi)')
% xlabel('\phi ')
% %print("Histograms_100_Bin_"+num2str(i),'-djpeg','-r500')
% %close all
% end
%  disp("Done")
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

function [HHSum,S1fTot] = getHistogramData(fname,Nel,Naz,sliceInds,pstruct)   
modTrue = pstruct.modTrue;
isLuna = pstruct.isLuna;
maxZ= pstruct.maxZ;
gpu = pstruct.gpu;
start = pstruct.start;
Nel = 1001;
Naz = 2001;
HHSum = zeros(Nel,Naz,9);
fractBins = 5;
fwx = 12;
dz = 5;
fwz = 1;
dopTh = 0.6;
wcorr = [];
rc = [];
dzres = 4.8;

for i = 1:size(sliceInds,2)

    sliceInd = sliceInds(i);
    disp("Working on Slice #"+num2str(i)+"/"+num2str(size(sliceInds,2))+"...")
    % Get S1
    %logF = replaceLogFile(fname,'    ','');
    st = struct;
    st.skipLastPoints = 30;
    st.window = 5;
    st.gpu = gpu;
    %st.logF = logF;

    if modTrue 
        if isLuna
            [t1,t2] = reconstructLunaWave(fname,sliceInd,fractBins);
            t(1,:,:,:) = squeeze(t1);
            t(2,:,:,:) = squeeze(t2);
        
            S1 = makeStokes(t(:,1:maxZ,:,:),3);
        else
            [S1,S2] = recstrTom(fname,sliceInd,st);

            if start > 1
                S1 = S2;
            end
        end
    else
        if isLuna
            [t1,t2] = reconstructLunaWave(fname,sliceInd,fractBins);
            t(1,:,:,:) = squeeze(t1);
            t(2,:,:,:) = squeeze(t2);
        else
            tom = recstrTom(fname,sliceInd,st);
            t(1,:,:,:) = squeeze(tom.ch1);
            t(2,:,:,:) = squeeze(tom.ch2);
        end
         S1 = makeStokes(t,3);
    end

    %% Filter & Get DOP
    dim = size(S1);

    if dim(4) == 3 % 3-component Stokes vector was provided
        S1 = cat(4,sqrt(dot(S1,S1,4)),S1);
    end 
    
    nx = (round(fwx*1.5)-1)/2;
    nx = linspace(-nx,nx,round(fwx*1.5))*2*sqrt(log(2))/fwx;
    h = exp(-nx.^2);
    if fwz>1
        nz = (round(fwaxial*1.5)-1)/2;
        nz = linspace(-nz,nz,round(fwaxial*1.5))*2*sqrt(log(2))/fwaxial;
        h = exp(-nz(:).^2)*h;
    end
    h = h/sum(h(:));
    S1f = imfilter(S1,h,'circular');
    %
    % final normalization
    QUVf = sqrt(dot(S1f(:,:,:,2:4),S1f(:,:,:,2:4),4));
    S1n = (S1f./QUVf);
    dop = mean(QUVf./S1f(:,:,:,1),3);
    mask = dop>0.9;
    % Remove Sheath
    mask(1:100,:) = 0;
    
    % Find Mirror State
    [HH,binsEl,binsAz] = poincareHistogram(S1n(:,:,:,2:4),Naz,Nel,mask);
    HHSum = HHSum+HH;
    S1fTot(:,:,:,:,i) = S1f;
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

% rowNum = floor(dim(1)/arraySize);
%     colNum = floor(dim(2)/arraySize);
%     % Cycle through pixels
%     for rowVal = 1:rowNum%dim(1)
% %         if(mod(row,10)==0)
% %             disp("Row = "+num2str(row))
% %         end
% 
%         for colVal = 1:colNum%dim(2)
%             % Determine if local maximum within 
%             row = (rowVal-1)*2*arraySize;
%             col = (colVal-1)*2*arraySize;
%             subArray = hBin(max(1,row-arraySize):min(dim(1),row+arraySize),max(1,col-arraySize):min(dim(2),col+arraySize));
%             % Determine maximum in this subarray
%             [maxVal] = max(max(subArray))
%             [minVal] = min(min(subArray));
%             if (maxVal-minVal) > minDiff
%                 % Record
%                 disp("Row - "+num2str(row)+", Col = "+num2str(col)+", N = "+num2str(nBin))
%                 disp("Diff = "+num2str(maxVal-minVal))
%                 [indx,indy]=find(subArray==maxVal)
% %                 figure,
% %                 subplot(1,2,1)
% %                 imagesc(hBin)
% %                 subplot(1,2,2)
% %                 imagesc(subArray)
% %                 hold on
% %                 plot(indy,indx,'sqr','MarkerSize',16,'LineWidth',2,'Color','r')
% %                 pause
%                 msMaxBin = [msMaxBin;row,col,maxVal-minVal];
%                 numMaxBin = numMaxBin+1;
%             end
%         end
%     end
% 
%     figure,
%     imagesc(hBin)
%     hold on
%     for i = 1:size(msMaxBin,1)
%        plot(msMaxBin(i,1),msMaxBin(i,2),'sqr','MarkerSize',16,'LineWidth',2,'Color','r') 
%        hold on
%     end
%     pause
% if numMaxBin>0
%     msMax{nBin} = msMaxBin;
%     numMax(nBin) = numMaxBin;
% else
%     maxVal = max(hBin(:));
%     [indx,indy]=find(subArray==maxVal)
%     msMax{nBin} = [indy,indx,0];
%     numMax(nBin) = 1;
% end

