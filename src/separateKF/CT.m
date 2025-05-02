clear;clc;close all

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

%% Constant Turn Rate Kalman Filter
w = deg2rad(5);
Act = [  eye(2) eye(2)*(sin(w*dt)/w) eye(2)*((1-cos(w*dt))/w);
       zeros(2)   eye(2)*(cos(w*dt))     eye(2)*(sin(w*dt)/w);
       zeros(2) eye(2)*(sin(w*dt)/w)        eye(2)*(cos(w*dt));];
B1 = [(sin(w*dt)/w) -((1-cos(w*dt))/w);((1-cos(w*dt))/w) (sin(w*dt)/w)];
B2 = [cos(w*dt) -sin(w*dt); sin(w*dt) cos(w*dt)];
Act = [  eye(2)       B1 zeros(2);
       zeros(2)       B2 zeros(2);
       zeros(2) zeros(2) zeros(2)];

qp = 1;
qv = 1;
qa = 10;
Qct = [qp*dt*eye(2)     zeros(2)     zeros(2);
           zeros(2) qv*dt*eye(2)     zeros(2);
           zeros(2)     zeros(2) qa*dt*eye(2);];


kfCT = kalmanFilter(Act,[],H,Qct,R,x0,P0);

for ii = 1 : n-1

    %% Run Kalman Filter 1
    % Time Update
    kfCT.propagate

    % Measurement Update
    kfCT.update(y(ii,:)')
    xCT_hat(ii,:) = kfCT.getState;

end

figure
plot(traj(:,1),traj(:,2))
hold on
plot(xCT_hat(:,1),xCT_hat(:,2))