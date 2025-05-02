clear;clc;close all

load("trajectoryData.mat")

n = length(traj);
xySigma  = 1;
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
type = "P";
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

%% Constant Velocity Kalman Filter
Acv = [  eye(2) dt*eye(2) zeros(2);
       zeros(2)    eye(2) zeros(2);
       zeros(2)  zeros(2) zeros(2);];
qp = 1;
qv = 1;
qa = 1;
Qcv = [qp*dt*eye(2)     zeros(2)     zeros(2);
           zeros(2) qv*dt*eye(2)     zeros(2);
           zeros(2)     zeros(2) qa*dt*eye(2);];

kfCV = kalmanFilter(Acv,[],H,Qcv,R,x0,P0);

for ii = 1 : n-1

    %% Run Kalman Filter 1
    % Time Update
    kfCV.propagate

    % Measurement Update
    kfCV.update(y(ii,:)')
    xCV_hat(ii,:) = kfCV.getState;

end

figure
plot(traj(:,1),traj(:,2))
hold on
plot(xCV_hat(:,1),xCV_hat(:,2))