function [LHS,RHS,outStruct] = findRotAndSheathCorrection(R1,R2)
%function [LHS,RHS] = findRotAndSheathCorrection(R1,R2) retrieves the LHS
%and RHS multipliers needed to be applied to the cumulative roundrtip 
%matrices to obtain absolute optic axis measurements. R1 and R1 are the
%cumulative round-trip signals from the inside and the outside of the
%birefringent sheath, which is assumed to be a Q-oriented linear retarder.
%
% We model the single pass transmission to the insight of the sheath as:
% R1fw = Rcath*Vth*L1, where L1 is a linear retarder, Vth is the
% deterministic rotation of the catheter, and Rcath is the retardance
% induced by the TIR of the probe and possibly additional shear from within
% the fiber.
%
% R2fw = Qsheath*R1fw, where Qsheath is the linear birefringence of the
% sheath, which varies as a function of rotation, even changing sign, but
% is assumed to be along the direction of Q.
%
% We first estimate L1, by 'centering' the sheath signal. Taking the root
% of the L1-corrected R1' provides:
% sqrt(R1') = Vth^-1*Lcath*Vth, where Lath is the 'linear' part of Rcath. 
%
% Computing sqrt(R1')^-1*R2'*sqrt(R1')^-1 then provides:
% (Vth^-1*Lcath^-1*Vth)*Vth^-1*D*RcathT*D*Qsheath*Qsheath*Rcath*Vth*(Vth^-1*Lcath^-1*Vth)
% = Vth^-1*Vcath^-1*Qsheath*Qsheath*Vcath*Vth
%
dim = size(R1);

% find L1
r1 = decomposeRot(R1);
fun = @(x) findL1(x,R1);
start = double(mean(r1,2)/2);
xopt = fminsearch(fun,start([1;2]));
[~,r1C] = fun(xopt);

Linv = makeRot3x3(cat(1,-xopt,0));
R1C = makeRot3x3(r1C);

% check for 'orientation' of ball lens signal
phi = unwrap(atan2(r1C(2,:),r1C(1,:)));
if mean(diff(phi))<0% negative rotation
    Linv = makeRot3x3(cat(1,(-xopt+pi*xopt/norm(xopt)),0));
    R1C = pagemtimes(Linv,pagemtimes(R1,Linv));
end


% just to make sure we have no wrap, unwrap this signal (in most instances
% this is unnecessary)
r1C = unwrapOA1D(R1C);
R1Csqinv = makeRot3x3(-r1C/2);

% correct R2 with the same Linv
R2C = pagemtimes(Linv,pagemtimes(R2,Linv));
%r2C = decomposeRot(R2C);

dRsheath = pagemtimes(R1Csqinv,pagemtimes(R2C,R1Csqinv));
drsheath = unwrapOA1D(dRsheath);

% also compute the 2pi wrapped version
drsheath2 = drsheath - 2*pi*drsheath./sqrt(sum(drsheath.^2,1));

%pick the one with the smaller ret
if mean(sqrt(sum(drsheath.^2,1)))>mean(sqrt(sum(drsheath2.^2,1)))
    drsheath = drsheath2;
end


Vtheta = makeRot3x3(cat(1,zeros(2,dim(3)),linspace(0,-4*pi/(dim(3)+1)*dim(3),dim(3))));
rloc = squeeze(pagemtimes(Vtheta,permute(drsheath,[1,3,2])));

% we want to find the theta offset that aligns the rotation vectors with
% the Q-axis. Considering SO3 matrices, this corresponds to matrices that
% are of the form [1,0,0;0,m22,m23;0,m32,m33]. Hence, we reduce the energy
% of the elements that we want to be zero:
fun = @(x)squeeze(sum(sum((makeRot(pagemtimes(makeRot3x3(cat(1,zeros(2,numel(x)),x)),permute(rloc,[1,3,4,2])))).^2.*[0;1;1;1;0;0;1;0;0],1),4));
angle1 = fminsearch(fun,0);
en1 = fun(angle1);

% this could be a local minimum; check optimim +pi/2
angle2 = fminsearch(fun,angle1+pi/2);
en2 = fun(angle2);%

if en2<en1
    meanAngle = angle2;
else
    meanAngle = angle1;
end

Vtheta = makeRot3x3(cat(1,zeros(2,dim(3)),linspace(0,-4*pi/(dim(3)+1)*dim(3),dim(3))+meanAngle));
rloc2 = squeeze(pagemtimes(Vtheta,permute(drsheath,[1,3,2])));

% make sure this is along positive Q
if mean(rloc2(1,:))<0
    meanAngle = meanAngle + pi;
    Vtheta = makeRot3x3(cat(1,zeros(2,dim(3)),linspace(0,-4*pi/(dim(3)+1)*dim(3),dim(3))+meanAngle));
%    rloc2 = squeeze(pagemtimes(Vtheta,permute(drsheath,[1,3,2])));
end

%fun = @(x)squeeze(sum(sum((makeRot(pagemtimes(makeRot3x3(cat(1,zeros(2,numel(x)),x)),permute(rloc,[1,2,4,3])))).^2.*[0;1;1;1;0;0;1;0;0],1),4));
%angleTemp = fminsearch(fun,0)

%
RHS = pagemtimes(Linv,pagemtimes(R1Csqinv,pagemtimes(makeRot3x3(-drsheath/2),pagetranspose(Vtheta))));
LHS = pagetranspose(RHS.*[1,1,-1;1,1,-1;-1,-1,1]);
%Tot = pagemtimes(pagetranspose(LHS),pagetranspose(RHS));

if nargout>2
    outStruct.meanAngle = meanAngle;
    outStruct.Linv = Linv;
    
    outStruct.r1 = decomposeRot(R1);
    outStruct.r2 = decomposeRot(R2);
    outStruct.r1C = decomposeRot(R1C);
    outStruct.r2C = decomposeRot(R2C);

    outStruct.drsheath = drsheath;

    outStruct.drsheathTheta = pagemtimes(Vtheta,permute(drsheath,[1,3,2]));
    outStruct.r1CTheta = pagemtimes(Vtheta,permute(r1C,[1,3,2]));
    outStruct.r2CTheta = pagemtimes(Vtheta,permute(outStruct.r2C,[1,3,2]));
end

%%
% 
% 
% Linv*R1Csqinv*Vth*qrsheath/2*Vth'
% 
% RHS = pagemtimes(makeRot3x3(drsheath/2),pagemtimes(makeRot3x3(r1C/2),pagemtimes((Vtheta),pagetranspose(Linv))));
% LHS = pagetranspose(RHS.*[1,1,-1;1,1,-1;-1,-1,1]);
% Tot = pagemtimes(LHS,RHS);
% 
% 
% 
% 
% 
% %R1C3 = pagemtimes(Linv,pagemtimes(R1,Linv));
% %r1C3 = decomposeRot(R1C3);
% 
% % check for 'orientation' of ball lens signal
% phi = unwrap(atan2(r1C(2,:),r1C(1,:)));
% if mean(diff(phi))<0% negative rotation
%     Linv = makeRot3x3(cat(1,(-xopt+2*pi*xopt/norm(xopt))/2));
%     R1C = pagemtimes(Linv,pagemtimes(R1,Linv));
%     r1C = decomposeRot(R1C);
% end
% % undo deterministic 4pi rotation
% Vtheta = makeRot3x3(cat(1,zeros(2,dim(3)),linspace(0,-4*pi/(dim(3)+1)*dim(3),dim(3))));
% 
% % recompute RsheathIn with Vtheta compensation
% r1C = squeeze(pagemtimes(Vtheta,permute(r1C,[1,4,2,3])));
% R1C = makeRot3x3(r1C);
% 
% R2C = pagemtimes(pagemtimes(Vtheta,Linv),pagemtimes(R2,pagemtimes(Linv,pagetranspose(Vtheta))));
% r2C = decomposeRot(R2C);
% 
% R1Csqinv = makeRot3x3(-r1C/2);
% 
% % get the retardance across the sheath
% dRsheath = pagemtimes(R1Csqinv,pagemtimes(R2C,R1Csqinv));
% drsheath = decomposeRot(dRsheath);
% drsheathuw = drsheath - 2*pi*drsheath./sqrt(sum(drsheath.^2,1));
% 
% % unwrap sheath; 
% mx = cat(2,zeros(3,1),decomposeRot(pagemtimes(makeRot3x3(-drsheath(:,1:end-1)/2),pagemtimes(makeRot3x3(drsheath(:,2:end)),makeRot3x3(-drsheath(:,1:end-1)/2)))));
% dd = cat(2,zeros(3,1),drsheath(:,2:end)-drsheath(:,1:end-1));
% 
% cp =(squeeze(sum(mx.*dd,1)./sqrt(sum(mx.^2,1))./sqrt(sum(dd.^2,1))));
% mm = mod(cumsum(cp<0),2)>0;
% %inds = cat(2,0,mod(cumsum(sum(mm(1:end-1)&~mm(2:end))>sum(mm(1:end-1)&mm(2:end),1)),2))>0;
% %mm(inds) = ~mm(inds);
% %mm = ~mm;
% 
% drsheath(:,~mm) = drsheathuw(:,~mm);
% 
% %Vtheta = squeeze(atan2(drsheath(2,:,:),drsheath(1,:,:)));
% % we want the signed non-resolved angle to align the sheath OA with the
% % Q-axis
% VV = makeRot3x3(cat(1,zeros(2,dim(3)),-atan(drsheath(2,:)./drsheath(1,:))));
% VVr = decomposeRot(VV);
% phiVV = squeeze(VVr(3,:));
% 
% Qsheath = pagemtimes(VV,pagemtimes(dRsheath,pagetranspose(VV)));
% qsheath = decomposeRot(Qsheath);
% 
% %%
% RHS = pagemtimes(makeRot3x3(drsheath/2),pagemtimes(makeRot3x3(r1C/2),pagemtimes((Vtheta),pagetranspose(Linv))));
% LHS = pagetranspose(RHS.*[1,1,-1;1,1,-1;-1,-1,1]);
% Tot = pagemtimes(LHS,RHS);
% 
% rtot = decomposeRot(Tot);
% 

visualize = false;
if visualize
%%
    r1 = decomposeRot(R1);
    r2 = decomposeRot(R2);

    r1C = decomposeRot(R1C);
    r2C = decomposeRot(R2C);

    drsheathTheta = pagemtimes(Vtheta,permute(drsheath,[1,3,2]));
    r1CTheta = pagemtimes(Vtheta,permute(r1C,[1,3,2]));
    r2CTheta = pagemtimes(Vtheta,permute(r2C,[1,3,2]));

    figure(2)
    clf
    subplot(2,2,1)
    plot(r1(1,:),r1(2,:),'.')
    hold on
    plot(r2(1,:),r2(2,:),'.')
    title('Rat sheath traces')
    
    subplot(2,2,2)
    plot(r1C(1,:),r1C(2,:),'.')
    hold on
    plot(r2C(1,:),r2C(2,:),'.')
    plot(drsheath(1,:),drsheath(2,:),'.')
    title('L1-compensated')
    legend('In','Out','Sheath')
    
    subplot(2,2,3)
    plot(r1CTheta(1,:),r1CTheta(2,:),'.')
    hold on
    plot(r2CTheta(1,:),r2CTheta(2,:),'.')
    plot(drsheathTheta(1,:),drsheathTheta(2,:),'.')
    title('Theta corrected')
    
    subplot(2,2,4)
    plot(atan(r1CTheta(2,:)./r1CTheta(1,:)),'.')
    hold on
    plot(atan(r2CTheta(2,:)./r2CTheta(1,:)),'.')
    plot(atan(drsheathTheta(2,:)./drsheathTheta(1,:)),'.')
%%
end


function [err,oacenter] = findL1(params,Mcath)
% [err,oacenter] = findL1(params,Mball) computes the error of the lsq trace
% minimization problem. oacenter is the optic axis of the catheter sheath
% signal after correcting for the linear element L1.

N = size(Mcath,3);

% Compensates Mcath with L1inv, which is a linear retarder
L1inv = makeRot3x3(-[params(1);params(2);0]);
Mcomp = pagemtimes(L1inv,pagemtimes(Mcath,L1inv));

% compute unwrapped ret and then compute the trace minimum
oacenter = decomposeRot(Mcomp); % unwrapping is necessary, as this is an artifact of the varying optic axis orientation; 
ss = repmat(cat(2,0,mod(cumsum(sqrt(sum(diff(oacenter,[],2).^2,1))>pi,2),2)),[3,1])>0;
delta = 2*pi*bsxfun(@rdivide,oacenter,sqrt(sum(oacenter.^2,1)));
oacenter(ss) = oacenter(ss) - delta(ss);
ret = sqrt(sum(oacenter.^2,1));
if mean(ret)>pi % the modulo 2*pi solutions do not change the sense of orientation, and simply switch the sign of mret
    oacenter = oacenter-2*pi*bsxfun(@rdivide,oacenter,sqrt(sum(oacenter.^2,1)));
    ret = sqrt(sum(oacenter.^2,1));
end
mret = atan2(mean(sin(real(ret))),mean(cos(real(ret))));
err = 1-mean(cos(ret-mret));

