clear;clc;close all

load("trajectoryData.mat")

n = length(traj);
xySigma  = 1;
VxySigma = 1;
axySigma = 1;
noise = [xySigma*randn([1 n]);
         VxySigma*randn([1 n]);
         axySigma*randn([1 n]);
         xySigma*randn([1 n]);
         VxySigma*randn([1 n]);
         axySigma*randn([1 n])];
meas = traj + noise';
y = meas(:,[1,4]);

tf = 60;
fs = 100;
dt = 1/fs;
time = 0 : dt : tf;
R = diag([VxySigma VxySigma]);

%% Kalman Filter 1
qCV = 400;
Q = qCV*dt*eye(3);
Q = [Q zeros(3);zeros(3) Q];
Q(3,3)=0;Q(6,6)=0;
A = [1  dt  0 0  0 0;
     0  1  0 0  0 0;
     0  0  0 0  0 0;
     0  0  0 1 dt 0;
     0  0  0 0  1 0;
     0  0  0 0  0 0];
H =  [1 0 0 0 0 0;
      0 0 0 1 0 0];
% H = eye(6);

x0 = [meas(1,:)]';
P0 = eye(6);

kf1 = kalmanFilter(A,[],H,Q,R,x0,P0);

%% Kalman Filter 2
qCT = 100;
Q = qCT*dt*eye(3);
Q = [Q zeros(3);zeros(3) Q];
Q(3,3)=0;Q(6,6)=0;

w = deg2rad(90/100);
A = [1    sin(w*dt)/w    0 0  -(1-cos(w*dt))/w 0;
      0      cos(w*dt)   0 0        -sin(w*dt) 0;
      0              0   0 0                 0 0
      0 (1-cos(w*dt))/w  0 1       sin(w*dt)/w 0;
      0       sin(w*dt)  0 0         cos(w*dt) 0;
      0               0  0 0                 0 0];

kf2 = kalmanFilter(A,[],H,Q,R,x0,P0);

%% Kalman Filter 3
qCA = 400;
Q = qCA*dt*eye(3);
Q = [Q zeros(3);zeros(3) Q];
A = [1  dt  .5*dt^2 0 0 0;
      0  1  dt 0 0 0
      0  0  1 0 0 0;
      0  0  0 1 dt .5*dt^2;
      0  0  0 0  1 dt;
      0  0  0 0  0 1];

kf3 = kalmanFilter(A,[],H,Q,R,x0,P0);

for ii = 1 : length(time)

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

figure
plot(traj(:,1),traj(:,4))
hold on
axis equal
plot(x1_hat(:,1),x1_hat(:,4),".")
plot(x2_hat(:,1),x2_hat(:,4),".")
plot(x3_hat(:,1),x3_hat(:,4),".")

figure 
plot(time,x1_hat(:,2))
hold on
plot(time,traj(:,2))
figure 
plot(time,x2_hat(:,2))
hold on
plot(time,traj(:,2))
figure 
plot(time,x3_hat(:,2))
hold on
plot(time,traj(:,2))

figure 
plot(time,x1_hat(:,5))
hold on
plot(time,traj(:,5))
figure 
plot(time,x2_hat(:,5))
hold on
plot(time,traj(:,5))
figure 
plot(time,x3_hat(:,5))
hold on
plot(time,traj(:,5))

figure 
plot(time,x1_hat(:,3))
hold on
plot(time,traj(:,3))
figure 
plot(time,x2_hat(:,3))
hold on
plot(time,traj(:,3))
figure 
plot(time,x3_hat(:,3))
hold on
plot(time,traj(:,3))

figure 
plot(time,x1_hat(:,6))
hold on
plot(time,traj(:,6))
figure 
plot(time,x2_hat(:,6))
hold on
plot(time,traj(:,6))
figure 
plot(time,x3_hat(:,6))
hold on
plot(time,traj(:,6))

% figure
% plot(time,x1_hat(:,2)-traj(:,3),".")
% figure
% plot(time,x2_hat(:,2)-traj(:,3),".")
% figure
% plot(time,x1_hat(:,2)-x2_hat(:,3),".")

% std(x1_hat(:,2)-traj(:,3))
% std(x2_hat(:,2)-traj(:,3))