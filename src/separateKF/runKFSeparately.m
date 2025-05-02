clear;clc;close all

load("trajectoryDataCVCTCA.mat")

n = length(traj);
xySigma  = .5;
VxySigma = .1;
axySigma = .01;
noise = [xySigma*randn([2 n]);
         VxySigma*randn([2 n]);
         axySigma*randn([2 n])];
meas = traj + noise';

fs = 100;
tf = (n-1)/fs;
dt = 1/fs;
time = 0 : dt : tf;

%% Kalman Filter Parameters
type = "PV";
switch type
    case "P"
        R = diag([xySigma xySigma])/dt;
        H = [1 0 0 0 0 0;
             0 1 0 0 0 0];
        y = meas(:,[1,2]);
        x0 = traj(1,:)';
        P0 = eye(6);
    case "V"
        R = diag([xySigma xySigma])/dt;
        H = [0 0 1 0 0 0;
             0 0 0 1 0 0];
        y = meas(:,[3,4]);
        x0 = traj(1,:)';
        P0 = eye(6);
    case "PV"
        R = diag([xySigma xySigma VxySigma VxySigma])/dt;
        H = [1 0 0 0 0 0;
             0 1 0 0 0 0;
             0 0 1 0 0 0;
             0 0 0 1 0 0];
        y = meas(:,[1:4]);
        x0 = traj(1,:)';
        P0 = eye(6);
    case "PVA"
        R = diag([xySigma xySigma VxySigma VxySigma AxySigma AxySigma])/dt;
        H = eye(6);
        y = meas;
        x0 = traj(1,:)';
        P0 = eye(6);
end


%% Kalman Filter 1
qCV = 10;
Q = qCV*dt*eye(3);
Q = [Q zeros(3);zeros(3) Q];
Q(3,3)=0;Q(6,6)=0;
A = [  eye(2) dt*eye(2) zeros(2);
       zeros(2)    eye(2) zeros(2);
       zeros(2)  zeros(2) zeros(2);];

kf1 = kalmanFilter(A,[],H,Q,R,x0,P0);

%% Kalman Filter 2
qCT = 10;
Q = qCT*dt*eye(3);
Q = [Q zeros(3);zeros(3) Q];
Q(3,3)=0;Q(6,6)=0;

w = deg2rad(30);
B1 = [(sin(w*dt)/w) -((1-cos(w*dt))/w);((1-cos(w*dt))/w) (sin(w*dt)/w)];
B2 = [cos(w*dt) -sin(w*dt); sin(w*dt) cos(w*dt)];
A = [  eye(2)       B1 zeros(2);
       zeros(2)       B2 zeros(2);
       zeros(2) zeros(2) zeros(2)];

kf2 = kalmanFilter(A,[],H,Q,R,x0,P0);

%% Kalman Filter 3
qCA = 1;
Q = qCA*dt*eye(3);
Q = [Q zeros(3);zeros(3) Q];
A = [  eye(2) dt*eye(2) .5*dt^2*eye(2);
       zeros(2)    eye(2) dt*eye(2);
       zeros(2)  zeros(2) eye(2);];

kf3 = kalmanFilter(A,[],H,Q,R,x0,P0);

for ii = 1 : length(time)-1

    %% Run Kalman Filter 1
    % Time Update
    kf1.propagate

    % Measurement Update
    kf1.update(y(ii,:)')
    x1_hat(ii,:) = kf1.getState;

    %% Run Kalman Filter 2
    % Time Update
    kf2.propagate

    % Measurement Update
    kf2.update(y(ii,:)')
    x2_hat(ii,:) = kf2.getState;

    %% Run Kalman Filter 3
    % Time Update
    kf3.propagate

    % Measurement Update
    kf3.update(y(ii,:)')
    x3_hat(ii,:) = kf3.getState;

end

fileName = 'D:\Classes\MECH 7710 - Optimal Control\Optimal-Final-Project\src\savedData\soloKF_CVCTCA';
save(fileName,'x1_hat','x2_hat','x3_hat')

figure
plot(x1_hat(:,3))
hold on
plot(x2_hat(:,3))
plot(x3_hat(:,3))
plot(traj(:,3))
% figure
% plot(time,x1_hat(:,2)-traj(:,3),".")
% figure
% plot(time,x2_hat(:,2)-traj(:,3),".")
% figure
% plot(time,x1_hat(:,2)-x2_hat(:,3),".")

% std(x1_hat(:,2)-traj(:,3))
% std(x2_hat(:,2)-traj(:,3))