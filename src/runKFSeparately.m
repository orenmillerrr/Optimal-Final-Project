clear;clc;close all

load("trajectory.mat")

n = length(X);
xySigma = 5;
VxySigma = .5;
noise = [xySigma*randn([1 n]);
         xySigma*randn([1 n]);
         VxySigma*randn([1 n]);
         VxySigma*randn([1 n])];
y = X + noise;

plot(X(1,:),X(3,:))
grid on
hold on
axis equal

tf = 200;
fs = 100;
dt = 1/fs;
time = 0 : dt : tf;
R = diag([xySigma xySigma VxySigma VxySigma]);

%% Kalman Filter 1

Q = eye(4);
A = [1 dt 0 0;
      0  1 0 0;
      0  0 1 dt;
      0  0 0  1];
H = [0 0 0 0;
      0 0 0 0;
      0 0 1 0;
      0 0 0 1];

x0 = X(:,1);
P0 = eye(4);

kf1 = kalmanFilter(A,[],H,Q,R,x0,P0);

%% Kalman Filter 2
w = 0.05;
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
    kf1.update(y(:,ii))

    x1_hat(:,ii) = kf1.getState;

    %% Run Kalman Filter 2
    % Time Update
    kf2.propagate

    % Measurement Update
    kf2.update(y(:,ii))

    x2_hat(:,ii) = kf2.getState;

end

plot(x1_hat(1,:),x1_hat(3,:),".")
plot(x2_hat(1,:),x2_hat(3,:),".")