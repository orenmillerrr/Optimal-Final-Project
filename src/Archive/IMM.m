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
Rd = diag([xySigma xySigma VxySigma VxySigma]);

%% Kalman Filter 1
P1 = eye(4);
Q1 = eye(4);
A1 = [1 dt 0 0;
      0  1 0 0;
      0  0 1 dt;
      0  0 0  1];
C1 = [0 0 0 0;
      0 0 0 0;
      0 0 1 0;
      0 0 0 1];
x_hat1 = X(:,1);

%% Kalman Filter 2
P2 = eye(4);
Q2 = eye(4);
w = 0.05;
A2 = [1    sin(w*dt)/w   0  -(1-cos(w*dt))/w;
      0      cos(w*dt)   0        -sin(w*dt);
      0 (1-cos(w*dt))/w   1      sin(w*dt)/w;
      0       sin(w*dt)  0         cos(w*dt)];
C2 = [0 0 0 0;
      0 0 0 0;
      0 0 1 0;
      0 0 0 1];
x_hat2 = X(:,1);

for ii = 1 : length(time)-1

    %% Run Kalman Filter 1
    % Time Update
    x_hat1_ = A1*x_hat1(:,ii);
    P1_ = A1*P1(:,:,ii)*A1' + Q1;

    % Measurement Update
    L1 = P1_*C1'*(C1*P1_*C1' + Rd)^-1;
    x_hat1(:,ii+1) = x_hat1_ + L1*(y(:,ii+1)-C1*x_hat1_);
    P1(:,:,ii+1) = (eye(size(P1_)) - L1*C1)*P1_;

    %% Run Kalman Filter 1
    % Time Update
    x_hat2_ = A2*x_hat2(:,ii);
    P2_ = A2*P2(:,:,ii)*A2' + Q2;

    % Measurement Update
    L2 = P2_*C2'*(C2*P2_*C2' + Rd)^-1;
    x_hat2(:,ii+1) = x_hat2_ + L2*(y(:,ii+1)-C2*x_hat2_);
    P2(:,:,ii+1) = (eye(size(P2_)) - L2*C2)*P2_;

end

plot(x_hat1(1,:),x_hat1(3,:),".")
plot(x_hat2(1,:),x_hat2(3,:),".")