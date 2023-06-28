%% Characterize System Compensation
function[ill_dgd,det_dgd,sps] = characterizeSysCom(sysComAv)
N = size(sysComAv.alignRotVec,2);
CSIPS = reshape(makeJones(sysComAv.symRotVec),[2,2,N]);
QTSIPS = reshape(makeJones(-sysComAv.alignRotVec),[2,2,N]);

% Q Alone
r = double(sysComAv.alignRotVec);
r2 =  double(sysComAv.symRotVec);
%CQT
for i = 1:N
   r2(:,i) = JonesDecomp(QTSIPS(:,:,i)*CSIPS(:,:,i));
end

N = size(sysComAv.alignRotVec,2);
Npoints = 1600;
deltaPhi = sqrt(sum((r(:,N)-r(:,1)).^2));%sqrt(sum(r(:,N).^2)) + sqrt(sum(r(:,1).^2))
deltaPhi2 = sqrt(sum((r2(:,N)-r2(:,1)).^2));%sqrt(sum(r(:,N).^2)) + sqrt(sum(r(:,1).^2))

qVec = r;
Minit = 0.5;
rlinInit= rand(3,1);

% Initial vector which will be optimized
vecInit = [rlinInit;Minit];
func = @(vec) lsqFunc(vec,qVec);
options = optimoptions(@lsqnonlin,'Display','none','FunctionTolerance',1e-12);
options.Algorithm = 'levenberg-marquardt';
[xFinal,err] = lsqnonlin(func,vecInit,[],[],options);
rlin = xFinal(1:3);
M = xFinal(4);
rVecFinal = (rlin*linspace(-M,M,N));%reshape(makeRot(rlin*linspace(0,1,N)),[3,3,N]);

lWidth = 1.5;
mSize = 20;
figure,
hold on
plot3(qVec(1,:),qVec(2,:),qVec(3,:),'r.-','MarkerSize',mSize,'LineWidth',lWidth)
plot3(rVecFinal(1,:),rVecFinal(2,:),rVecFinal(3,:),'b.-','MarkerSize',mSize,'LineWidth',lWidth)
[x,y,z] = sphere;
mesh(x,y,z,'FaceColor','none','EdgeColor',[0.3,0.3,0.3])
xlabel('Q','FontSize',20,'Color','k')
ylabel('U','FontSize',20,'Color','k')
zlabel('V','FontSize',20,'Color','k')
set(gca,'xtick',[-1,1],'ytick',[-1,1],'ztick',[-1,1]);
legend("Ilumination System Compensation","Model Fit")
deltaPhiQAlt = (M*2)*1.2;

cVec = r2;
Minit = 0.5;
rlinInit= rand(3,1);
rdcInit= rand(3,1);

% Initial vector which will be optimized
vecInit = [rlinInit;rdcInit;Minit];
func = @(vec) lsqFunc2(vec,cVec);
options = optimoptions(@lsqnonlin,'Display','none','FunctionTolerance',1e-12);
options.Algorithm = 'levenberg-marquardt';
[xFinal,err] = lsqnonlin(func,vecInit,[],[],options);
rdcFinal = xFinal(1:3);
rlinFinal = xFinal(4:6);
Mfinal = xFinal(7);
RDC = makeRot(rdcFinal);%reshape(makeRot(rdc),[3,3]);
Rlin = makeRot(rlinFinal*linspace(-Mfinal,Mfinal,N));%reshape(makeRot(rlin*linspace(0,1,N)),[3,3,N]);
rVecFinal = decomposeRot(MatrixMultiply(RDC,Rlin));%pagemtimes(RDC,Rlin);
    
lWidth = 1.5;
mSize = 20;

figure,
hold on
plot3(cVec(1,:),cVec(2,:),cVec(3,:),'r.-','MarkerSize',mSize,'LineWidth',lWidth)
plot3(rVecFinal(1,:),rVecFinal(2,:),rVecFinal(3,:),'b.-','MarkerSize',mSize,'LineWidth',lWidth)
[x,y,z] = sphere;
mesh(x,y,z,'FaceColor','none','EdgeColor',[0.3,0.3,0.3])
xlabel('Q','FontSize',20,'Color','k')
ylabel('U','FontSize',20,'Color','k')
zlabel('V','FontSize',20,'Color','k')
set(gca,'xtick',[-1,1],'ytick',[-1,1],'ztick',[-1,1])
legend("Detection System Compensation","Model Fit")

% Find Illumination DGD
%Add 20% onto spectrum
Nk = Npoints/2;
Nfft = 32768;
Nzeropad = Nfft-Nk;
scale = Nk/Nfft;

% Find illumination DGD
e1 = exp(1i*(0:Nk-1)*deltaPhi/Nk);
e1_zp = cat(1,e1.',zeros(Nzeropad,1));
f1 = fft(e1_zp);
[peaks,peaksLoc] = findpeaks(abs(f1));
[maxVal,maxLoc] = max(peaks);
maxPeakLoc = peaksLoc(maxLoc);
ill_dgd = 5.61*scale*maxPeakLoc;

% Find detection DGD
e1 = exp(1i*(0:Nk-1)*deltaPhi2/Nk);
e1_zp = cat(1,e1.',zeros(Nzeropad,1));
f1 = fft(e1_zp);
[peaks,peaksLoc] = findpeaks(abs(f1));
[maxVal,maxLoc] = max(peaks);
maxPeakLoc = peaksLoc(maxLoc);
det_dgd = 5.61*scale*maxPeakLoc;

% Find SPS
Mq = makeRot3x3(r);
bp = squeeze(Mq(1,:,:));
diffbp = bp(:,2:end)-bp(:,1:end-1);
sps = sum(sqrt(sum(diffbp.^2,1)));
end
function[err] = lsqFunc(vec,qVec)
N = 9;
rlin = vec(1:3);
M = vec(4);
rVec = (rlin*linspace(-M,M,N));%reshape(makeRot(rlin*linspace(0,1,N)),[3,3,N]);
 
% Compute LSQ
 err = 0;
 for i = 1:N
    binErr = sum(abs(qVec(:,i)-rVec(:,i)).^2);
    err = err + binErr;
 end

 
end

function[err] = lsqFunc2(vec,points)
N = 9;
rdc = vec(1:3);
rlin = vec(4:6);
 M = vec(7);
 RDC = makeRot(rdc);%reshape(makeRot(rdc),[3,3]);
 Rlin = makeRot(rlin*linspace(-M,M,N));%reshape(makeRot(rlin*linspace(0,1,N)),[3,3,N]);
 rVec = decomposeRot(MatrixMultiply(RDC,Rlin));%pagemtimes(RDC,Rlin);
 

 % Compute LSQ
 err = 0;
 for i = 1:N
    binErr = sum(abs(points(:,i)-rVec(:,i)).^2);
    err = err + binErr;
 end

 
end
