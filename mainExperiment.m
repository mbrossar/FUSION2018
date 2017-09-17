%data comes from the EUROC MAV dataset V1_01_medium
clear;
close all;
addpath myToolbox
addpath filters
addpath data

tMin = 1081; % starting IMU time
tImagesMin = 109; % starting camera time

%% Initialization
freqIMU = 200; %Hz
freqCam = 20; %Hz

ParamGlobal.dirImage = 'V1_02_medium/mav0/cam0/data/';
dirImage = 'V1_02_medium/mav0/cam0/data/';
fileData = 'DATA_MEDIUM.mat';
fileImages = 'fileImages_MEDIUM.mat';
offset = 881;

% IMU and frame date
load(fileData);
t = tIMU;
t = t/10^9; % ns -> s
g = 9.81*[0;0;-1]; %gravity field
omega = omega(:,tMin:end);
acc = acc(:,tMin:end);
t = t(tMin:end);
tIMU = tIMU(tMin:end);
ParamGlobal.tIMU = tIMU;
NbSteps = length(tIMU);
NbStepsMax = 8000;

load(fileImages); 
tImages = tImages(tImagesMin:end);
fileImages = fileImages(tImagesMin:end);
ParamGlobal.fileImages = fileImages;
%corect offset
tReal(1:offset) = [];
trajReal.x(:,1:offset) = []; 
trajReal.v(:,1:offset) = [];
trajReal.quat(1:offset) = [];
trajReal.psi(1:offset) = [];
trajReal.theta(1:offset) = [];
trajReal.phi(1:offset) = [];
trajReal.omega_b(:,1:offset) = [];
trajReal.a_b(:,1:offset) = [];
obsTimes = zeros(length(t),1);
obsTimes(1:10:end) = 1; % IMU is 10 times faster than camera

%% Init Filter
P0amers = diag([1;1;1]*1.e-3); %initial landmark covariance
R = 1.0^2*eye(2); %measurement noise for one landmark
ParamFilter.NbAmers = 30; %nominal number of landmarks in the state
ParamFilter.NbAmersMin = ParamFilter.NbAmers;
ParamFilter.EcartPixelMax = 20;

% init covariance
p0Rot = (0.01*pi/180)^2;
p0v =  1.e-4;
p0x =  1.e-8;
p0omegab = 1.e-6;
p0ab = 1.e-6;
P0 = diag([p0Rot*ones(3,1);p0v*ones(3,1);p0x*ones(3,1);...
    p0omegab*ones(3,1);p0ab*ones(3,1)]);

% process noises
q_omega = (1.6968e-4)^2*200;
q_a = (2e-3)^2*200;
q_omegab = (1.9393e-5)^2*200;
q_ab = (3e-3)^2*200;
Q = diag([q_omega*ones(3,1);q_a*ones(3,1);q_omegab* ...
    ones(3,1);q_ab*ones(3,1)]);
Qc = chol(Q);
P0 = blkdiag(P0,kron(eye(ParamFilter.NbAmers),P0amers));

%depending on the chosen camera %cam0
ParamFilter.chiC = [0.0148655429818, -0.999880929698, 0.00414029679422, -0.0216401454975;
    0.999557249008, 0.0149672133247, 0.025715529948, -0.064676986768;
    -0.0257744366974, 0.00375618835797, 0.999660727178, 0.00981073058949;
    0.0, 0.0, 0.0, 1.0]; % camera pose
ParamFilter.Pi = [458.654 0 0; 0 457.296 0; 367.215, 248.375 1]';%camera calibration matrix
ParamFilter.cameraParams = cameraParameters('IntrinsicMatrix', ParamFilter.Pi',...
    'RadialDistortion',[-0.28340811, 0.07395907],...
    'TangentialDistortion',[0.00019359, 1.76187114e-05]);

% Initialisation of the state is obtained following [Mur-Artal,2017],
load('data/ORB_SLAM_init.mat');
Rot0 = eul2rotm([trajReal.psi(1),trajReal.theta(1),trajReal.phi(1)]);
x0 = trajReal.x(:,1);
v0 = trajReal.v(:,1);
omega_b0 = orb_slam.omega_b;
a_b0 = orb_slam.a_b;
PosAmers0 = orb_slam.PosAmers;
trackerMain = orb_slam.trackerMain;
trackerBis = orb_slam.trackerBis;
myTracks = orb_slam.myTracks;

IdxImage = 2; % image index

trajL = initTraj(NbSteps);
RotL = Rot0;
vL = v0;
xL = x0;
omega_bL = omega_b0;
a_bL = a_b0;
PosAmersL = PosAmers0;
P_L = P0;
S_L = chol(P_L);
chiL = state2chi(RotL,vL,xL,PosAmersL);

trajR = initTraj(NbSteps);
RotR = Rot0;
vR = v0;
xR = x0;
omega_bR = omega_b0;
a_bR = a_b0;
PosAmersR = PosAmers0;
P_R = P0;
S_R = chol(P_R);
chiR = state2chi(RotR,vR,xR,PosAmersR);

trajRef = initTraj(NbSteps);
RotRef = Rot0;
vRef = v0;
xRef = x0;
omega_bRef = omega_b0;
a_bRef = a_b0;
PosAmersRef = PosAmers0;
P_Ref = blkdiag(P0);
S_Ref = chol(P_Ref);
chiRef = [RotRef xRef;0 0 0 1];
xidotRef = zeros(6,1);

trajU = initTraj(NbSteps);
RotU = Rot0;
vU = v0;
xU = x0;
omega_bU = omega_b0;
a_bU = a_b0;
PosAmersU = PosAmers0;
P_U = blkdiag(P0);
S_U = chol(P_U);

trajI = initTraj(NbSteps);
RotI = Rot0;
vI = v0;
xI = x0;
omega_bI = omega_b0;
a_bI = a_b0;
PosAmersI = PosAmers0;
P_I = blkdiag(P0);

%% Filtering
for i = 2:NbStepsMax
    % propagation
    dt = t(i)-t(i-1);
    omega_i = omega(:,i);
    acc_i = acc(:,i);
    
    chiAntR = chiR;
    RotR = RotR*expSO3((omega_i-omega_bR)*dt);
    vR = vR+(RotR*(acc_i-a_bR)+g)*dt;
    xR = xR+vR*dt;
    chiR = state2chi(RotR,vR,xR,PosAmersR);
    S_R = rukfPropagation(dt,chiR,chiAntR,omega_bR,a_bR,S_R,omega_i,...
        acc_i,Qc,g);
    
    % propagation for others filters
    chiAntL = chiL;
    RotL = RotL*expSO3((omega_i-omega_bL)*dt);
    vL = vL+(RotL*(acc_i-a_bL)+g)*dt;
    xL = xL+vL*dt;
    chiL = state2chi(RotL,vL,xL,PosAmersL);
    S_L = lukfPropagation(dt,chiL,chiAntL,omega_bL,a_bL,S_L,omega_i,...
        acc_i,Qc,g);
    
    chiAntRef = [RotRef,xRef;zeros(1,3) 1];
    vAnt = vRef;
    RotRef = RotRef*expSO3((omega_i-omega_bRef)*dt);
    vRef = vRef+(RotRef*(acc_i-a_bRef)+g)*dt;
    xRef = xRef+vRef*dt;
    xidotRef = xidotRef + [omega_i-omega_bRef;vRef]*dt;
    chiRef = [RotRef,xRef;zeros(1,3) 1];
    S_Ref = ukfRefPropagation(dt,chiRef,chiAntRef,vAnt,omega_bRef,...
        a_bRef,S_Ref,omega_i,acc_i,Qc,g);
    
    RotAntU = RotU;
    RotU = RotU*expSO3((omega_i-omega_bU)*dt);
    vU = vU+(RotU*(acc_i-a_bU)+g)*dt;
    xU = xU+vU*dt;
    chiU = [RotU,xU;zeros(1,3) 1];
    S_U = ukfPropagation(dt,RotU,RotAntU,vAnt,omega_bU,a_bU,S_U,omega_i,...
        acc_i,Qc,g);
    
    [RotI,vI,xI,PosAmersI,P_I] = iekfPropagation(dt,RotI,vI,xI,...
        omega_bI,a_bI,PosAmersI,P_I,omega_i,acc_i,Q,g);
    chiI = state2chi(RotI,vI,xI,PosAmersI);
    
    % if measurement
    if obsTimes(i) == 1
        % track points in image
        [y,yAmers,trackerMain,trackerBis,pointsMain,validityMain,...
            myTracks,pointsBis] = ...
            ObserveLandmarks(trackerMain,trackerBis,dirImage,IdxImage,...
            fileImages,ParamFilter,RotR,xR,PosAmersR,i,S_R,myTracks);
        
        % update state
        param.yAmers = yAmers;
        [chiR,omega_bR,a_bR,S_R] = rukfUpdate(chiR,omega_bR,...
            a_bR,S_R,y,param,R,ParamFilter);
        [RotR,vR,xR,PosAmersR] = chi2state(chiR);
        
        % update other filters
        chiL = state2chi(RotL,vL,xL,PosAmersL);
        [chiL,omega_bL,a_bL,S_L] = lukfUpdate(chiL,omega_bL,...
            a_bL,S_L,y,param,R,ParamFilter);
        [RotL,vL,xL,PosAmersL] = chi2state(chiL);
        
        param.PosAmers = PosAmersRef;
        chiRef = state2chi(RotRef,vRef,xRef,PosAmersRef);
        [chiRef,vRef,PosAmersRef,omega_bRef,a_bRef,S_Ref,xidotRef] = ukfRefUpdate(chiRef,vRef,omega_bRef,...
            a_bRef,S_Ref,y,param,R,ParamFilter,PosAmersRef,xidotRef);
        RotRef = chiRef(1:3,1:3);
        xRef = chiRef(1:3,4);
        
        param.PosAmers = PosAmersU;
        chiU = state2chi(RotU,vU,xU,PosAmersU);
        [RotU,vU,xU,PosAmersU,omega_bU,a_bU,S_U] = ukfUpdate(RotU,vU,xU,omega_bU,...
            a_bU,S_U,y,param,R,ParamFilter,PosAmersU);
        
        [chiI,omega_bI,a_bI,P_I] = iekfUpdate(chiI,omega_bI,a_bI,...
            P_I,y,R,ParamFilter,yAmers);
        [RotI,vI,xI,PosAmersI] = chi2state(chiI);
        
        % save trajectory
        trajR = updateTraj(trajR,RotR,vR,xR,omega_bR,a_bR,i);
        trajL = updateTraj(trajL,RotL,vL,xL,omega_bL,a_bL,i);
        trajRef = updateTraj(trajRef,RotRef,vRef,xRef,omega_bRef,a_bRef,i);
        trajU = updateTraj(trajU,RotU,vU,xU,omega_bU,a_bU,i);
        trajI = updateTraj(trajI,RotI,vI,xI,omega_bI,a_bI,i);
        
        % remplace non visible landmarks
        [S_R,PosAmersR,ParamFilter,trackerBis,myTracks,PosAmersNew,...
            IdxAmersNew,trackCov,pointsMain,validityMain] = manageAmers(S_R,...
            PosAmersR,ParamFilter,ParamGlobal,trackerBis,...
            trajR,i,pointsMain,validityMain,IdxImage,myTracks,pointsBis);
        chiR = state2chi(RotR,vR,xR,PosAmersR);
        
        setPoints(trackerMain,pointsMain,validityMain);
        
        % remplace non visible landmarks for others filters with same
        % landmarks
        if isempty(IdxAmersNew) == 0
            P_L = S_L'*S_L;
            P_Ref = S_Ref'*S_Ref;
            P_U = S_U'*S_U;
            for jj = 1:length(IdxAmersNew)
                idx = IdxAmersNew(jj);
                idxP = 15+(3*idx-2:3*idx);
                P_L(:,idxP) = 0;
                P_L(idxP,:) = 0;
                P_L(idxP,idxP) = trackCov{jj};
                PosAmersL(:,idx) = PosAmersNew(jj,:)';
                P_Ref(:,idxP) = 0;
                P_Ref(idxP,:) = 0;
                P_Ref(idxP,idxP) = trackCov{jj};
                PosAmersRef(:,idx) = PosAmersNew(jj,:)';
                P_U(:,idxP) = 0;
                P_U(idxP,:) = 0;
                P_U(idxP,idxP) = trackCov{jj};
                PosAmersU(:,idx) = PosAmersNew(jj,:)';
                P_I(:,idxP) = 0;
                P_I(idxP,:) = 0;
                P_I(idxP,idxP) = trackCov{jj};
                PosAmersI(:,idx) = PosAmersNew(jj,:)';
            end
            S_L = chol(P_L);
            chiL = state2chi(RotL,vL,xL,PosAmersL);
            S_Ref = chol(P_Ref);
            chiRef = [RotRef,xRef;zeros(1,3) 1];
            S_U = chol(P_U);
            chiI = state2chi(RotI,vI,xI,PosAmersI);
        end
        disp(i/200);
        IdxImage = IdxImage+1;
    else
        trajR = updateTraj(trajR,RotR,vR,xR,omega_bR,a_bR,i);
        trajL = updateTraj(trajL,RotL,vL,xL,omega_bL,a_bL,i);
        trajRef = updateTraj(trajRef,RotRef,vRef,xRef,omega_bRef,a_bRef,i);
        trajU = updateTraj(trajU,RotU,vU,xU,omega_bU,a_bU,i);
        trajI = updateTraj(trajI,RotI,vI,xI,omega_bI,a_bI,i);
    end
end

%% Plots
errorR = computeError(trajR,trajReal,i);
errorL = computeError(trajL,trajReal,i);
errorU = computeError(trajU,trajReal,i);
errorRef = computeError(trajRef,trajReal,i);
errorI = computeError(trajI,trajReal,i);

figure;hold on;
plot(t(2:i-1)-t(2),errorR.errorR);
plot(t(2:i-1)-t(2),errorL.errorR);
plot(t(2:i-1)-t(2),errorRef.errorR);
plot(t(2:i-1)-t(2),errorU.errorR);
plot(t(2:i-1)-t(2),errorI.errorR);
disp(sqrt(mean([errorR.errorR errorL.errorR  errorRef.errorR errorU.errorR errorI.errorR].^2)));
legend('R-UKF-LG','L-UKF-LG','SE(3)-UKF','UKF','IEKF')
xlabel('t (s)')
ylabel('RMSE attitdude (°)')
title('RMSE on attitude as function of time')
figure
hold on;
plot(t(2:i-1)-t(2),errorR.errorX);
plot(t(2:i-1)-t(2),errorL.errorX);
plot(t(2:i-1)-t(2),errorRef.errorX);
plot(t(2:i-1)-t(2),errorU.errorX);
plot(t(2:i-1)-t(2),errorI.errorX);
disp(sqrt(mean([errorR.errorX  errorL.errorX errorRef.errorX errorU.errorX errorI.errorX].^2)));
legend('R-UKF-LG','L-UKF-LG','SE(3)-UKF','UKF','IEKF')
xlabel('t (s)')
ylabel('RMSE position (m)')
title('RMSE on position as function of time')