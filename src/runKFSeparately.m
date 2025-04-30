clear;clc;close all

load("trajectory_data.mat")

n = length(traj);
xySigma = 1;
VxySigma = 1;
axySigma = .1;
noise = [xySigma*randn([2 n]);
         VxySigma*randn([2 n]);
         axySigma*randn([2 n]);];
y = traj + noise';
y = y(:,[3,4]);

plot(traj(:,1),traj(:,2),".")
grid on
hold on
axis equal

tf = 65;
fs = 100;
dt = 1/fs;
time = 0 : dt : tf;
R = diag([VxySigma VxySigma]);

%% Kalman Filter 1
qCV = 400;
Q = qCV*dt*[dt^2/3 dt/2;dt/2 1];
Q = [Q zeros(2);zeros(2) Q];
A = [1 dt 0 0;
      0  1 0 0;
      0  0 1 dt;
      0  0 0  1];
H =  [0 1 0 0;
      0 0 0 1];

x0 = [traj(1,1) traj(3,1) traj(2,1) traj(4,1)]';
P0 = eye(4);

kf1 = kalmanFilter(A,[],H,Q,R,x0,P0);

%% Kalman Filter 2
qTURN = 25;
Q = qCV*dt*[1 0;0 1];
Q = [Q zeros(2);zeros(2) Q];

w = deg2rad(30);
A = [1    sin(w*dt)/w   0  -(1-cos(w*dt))/w;
      0      cos(w*dt)   0        -sin(w*dt);
      0 (1-cos(w*dt))/w   1      sin(w*dt)/w;
      0       sin(w*dt)  0         cos(w*dt)];

kf2 = kalmanFilter(A,[],H,Q,R,x0,P0);

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

end

plot(x1_hat(:,1),x1_hat(:,3),".")
plot(x2_hat(:,1),x2_hat(:,3),".")

% figure
% plot(time,x1_hat(:,2)-traj(:,3),".")
% figure
% plot(time,x2_hat(:,2)-traj(:,3),".")
% figure
% plot(time,x1_hat(:,2)-x2_hat(:,3),".")

% std(x1_hat(:,2)-traj(:,3))
% std(x2_hat(:,2)-traj(:,3))