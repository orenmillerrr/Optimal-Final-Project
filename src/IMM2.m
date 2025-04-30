clear;clc;close all

load("trajectory/trajectoryData.mat")

n = length(traj);
xySigma = 1;
VxySigma = 1;
axySigma = .1;
noise = [xySigma*randn([1 n]);
         VxySigma*randn([1 n]);
         axySigma*randn([1 n]);
         xySigma*randn([1 n]);
         VxySigma*randn([1 n]);
         axySigma*randn([1 n])];
y = traj + noise';
y = y(:,[2,5]);

tf = 60;
fs = 100;
dt = 1/fs;
time = 0 : dt : tf;
R = diag([VxySigma VxySigma]);

% %% Kalman Filter 1
% qCV = 400;
% Q = qCV*dt*[dt^2/3 dt/2;dt/2 1];
% Q = [Q zeros(2);zeros(2) Q];
% A = [1 dt 0 0;
%      0  1 0 0;
%      0  0 1 dt;
%      0  0 0  1];
% H = [0 1 0 0;
%      0 0 0 1];
% 
% x0 = [traj(1,1) traj(1,2) traj(1,3) traj(1,4)]';
% P0 = eye(4);
% 
% kf1 = kalmanFilter(A,[],H,Q,R,x0,P0);

%% Kalman Filter 1
qCV = 400;
Q = qCV*dt*[dt^2/3 dt/2 0;dt/2 1 0;0 0 0];
Q = [Q,zeros(3);zeros(3),Q];
A = [1 dt 0 0  0 0;
     0  1 0 0  0 0;
     0  0 0 0  0 0;
     0  0 0 1 dt 0;
     0  0 0 0  1 0;
     0  0 0 0  0 0];

H = [0 1 0 0 0 0;
     0 0 0 0 1 0];

x0 = traj(1,:)';
P0 = eye(6);

kf1 = kalmanFilter(A,[],H,Q,R,x0,P0);

%% Kalman Filter 2
qTURN = 25;
Q = qCV*dt*eye(3);
Q = [Q zeros(3);zeros(3) Q];

w = deg2rad(90/5);
A = [1     sin(w*dt)/w   (1-cos(w*dt))/(w^2) 0 0 0;
     0       cos(w*dt)    sin(w*dt)/w           0 0 0;
     0      -w*sin(w*dt) cos(w*dt)              0 0 0;
     0 0 0 1     sin(w*dt)/w   (1-cos(w*dt))/(w^2);
     0 0 0 0     cos(w*dt)     sin(w*dt)/w;
     0 0 0 0     -w*sin(w*dt) cos(w*dt)];

kf2 = kalmanFilter(A,[],H,Q,R,x0,P0);

%% Kalman Filter 3
qCA = 400;
Q = qCA*dt*[0 0 0;0 0 0;0 0 1];
Q = [Q,zeros(3);zeros(3),Q];

A = [1 dt .5*dt^2 0  0 0;
     0  1 dt 0  0 0;
     0  0 1 0  0 0;
     0  0 0 1 dt .5*dt^2;
     0  0 0 0  1 dt;
     0  0 0 0  0 1];


kf3 = kalmanFilter(A,[],H,Q,R,x0,P0);

%% IMM Parameters
x{1}  = x0;
P{1}  = P0;
x{2}  = x0;
P{2}  = P0;
x{3}  = x0;
P{3}  = P0;
M = [0.90 0.05 0.05;
      0.05 0.90 0.05;    % Markov Chain
      0.05 0.05 0.10];    
mu = [1/3 1/3 1/3];
cBar = mu * M;
numModels = 3;

for ii = 1 : length(time)

    [Z1,S1] = kf1.computeInnovation(y(ii,:)');
    [L1] = kf1.computeLikelihood(Z1,S1);

    [Z2,S2] = kf2.computeInnovation(y(ii,:)');
    [L2] = kf2.computeLikelihood(Z2,S2);

    [Z3,S3] = kf3.computeInnovation(y(ii,:)');
    [L3] = kf3.computeLikelihood(Z3,S3);

    % IMM Likelyhood
    likeIMM = [L1 L2 L3];

    % IMM Weights
    weightsIMM = cBar(ii,:) .* likeIMM;
    weightsIMM = weightsIMM/sum(weightsIMM);

    % Measurement Update
    kf1.update(y(ii,:)')
    x1_hat(:,ii) = kf1.getState;
    x{1} = x1_hat(:,ii);
    P{1} = kf1.getCovariance;

    % Measurement Update
    kf2.update(y(ii,:)')
    x2_hat(:,ii) = kf2.getState;
    x{2} = x2_hat(:,ii);
    P{2} = kf2.getCovariance;

    % Measurement Update
    kf3.update(y(ii,:)')
    x3_hat(:,ii) = kf3.getState;
    x{3} = x3_hat(:,ii);
    P{3} = kf3.getCovariance;

    xIMM(ii,:) = weightsIMM(1)*x{1} + weightsIMM(2)*x{2} + weightsIMM(3)*x{3};

    pIMM = (P{1} + (x{1} -  xIMM(ii,:))'*(x{1} -  xIMM(ii,:)))*weightsIMM(1) + ...
              (P{2} + (x{2} -  xIMM(ii,:))'*(x{2} -  xIMM(ii,:)))*weightsIMM(2) + ...
              (P{3} + (x{3} -  xIMM(ii,:))'*(x{3} -  xIMM(ii,:)))*weightsIMM(3);

    % Compute Mixed Initial Conditions
    cBar1= M(:,1)'*weightsIMM';  
    if cBar1>1e-80
        omega(1,1)=M(1,1)*weightsIMM(1)/cBar1;
        omega(2,1)=M(2,1)*weightsIMM(2)/cBar1;
        omega(3,1)=M(3,1)*weightsIMM(3)/cBar1;
    else
        cBar1 = 0;
        omega(1,1) = 0;
        omega(2,1) = 0;
        omega(3,1) = 0;
    end
    cBar2=M(:,2)'*weightsIMM';
    if cBar2>1e-80
        omega(1,2)=M(1,2)*weightsIMM(1)/cBar2;
        omega(2,2)=M(2,2)*weightsIMM(2)/cBar2;
        omega(3,2)=M(3,2)*weightsIMM(3)/cBar2;
    else
        cBar2 = 0;
        omega(1,2) = 0;
        omega(2,2) = 0;
        omega(3,2) = 0;
    end
    cBar3=M(:,3)'*weightsIMM';
    if cBar3>1e-80
        omega(1,3)=M(1,3)*weightsIMM(1)/cBar3;
        omega(2,3)=M(2,3)*weightsIMM(2)/cBar3;
        omega(3,3)=M(3,3)*weightsIMM(3)/cBar3;
    else
        cBar3 = 0;
        omega(1,3) = 0;
        omega(2,3) = 0;
        omega(3,3) = 0;
    end

    cBar(ii+1,:)=[cBar1,cBar2, cBar3]; 

    x1Mixed = omega(1,1)*x{1}+omega(2,1)*x{2}+omega(3,1)*x{3};
    x2Mixed = omega(1,2)*x{1}+omega(2,2)*x{2}+omega(3,2)*x{3};
    x3Mixed = omega(1,3)*x{1}+omega(2,3)*x{2}+omega(3,3)*x{3};

    P1Mixed = (kf1.P + (kf1.x_hat - x1Mixed)'*(kf1.x_hat - x1Mixed))*omega(1,1) + ...
                 (kf2.P + (kf2.x_hat - x2Mixed)'*(kf2.x_hat - x2Mixed))*omega(2,1) + ...
                 (kf3.P + (kf3.x_hat - x3Mixed)'*(kf3.x_hat - x3Mixed))*omega(3,1); 

    P2Mixed = (kf1.P + (kf1.x_hat - x1Mixed)'*(kf1.x_hat - x1Mixed))*omega(1,2) + ...
                 (kf2.P + (kf2.x_hat - x2Mixed)'*(kf2.x_hat - x2Mixed))*omega(2,2) + ...
                 (kf3.P + (kf3.x_hat - x3Mixed)'*(kf3.x_hat - x3Mixed))*omega(3,2);

    P3Mixed = (kf1.P + (kf1.x_hat - x1Mixed)'*(kf1.x_hat - x1Mixed))*omega(1,3) + ...
                 (kf2.P + (kf2.x_hat - x2Mixed)'*(kf2.x_hat - x2Mixed))*omega(2,3) + ...
                 (kf3.P + (kf3.x_hat - x3Mixed)'*(kf3.x_hat - x3Mixed))*omega(3,3);

    kf1.x_hat = x1Mixed;
    kf1.P = P1Mixed;

    kf2.x_hat = x2Mixed;
    kf2.P = P2Mixed;

    kf3.x_hat = x3Mixed;
    kf3.P = P3Mixed;

    % Time Update
    kf1.propagate
    kf2.propagate
    kf3.propagate
end

%% Plot Trajectory Figures
% XY Positions
plot(time,traj(:,1),".")
grid on
hold on
plot(time,traj(:,4),".")
plot(time,vecnorm([traj(:,1),traj(:,4)]'),"--k")

% XY Velocities
figure
plot(time,traj(:,2),".")
grid on
hold on
plot(time,traj(:,5),".")
plot(time,vecnorm([traj(:,2),traj(:,5)]'),"--k")

% XY Accelerations
figure
plot(time,traj(:,3),".")
grid on
hold on
plot(time,traj(:,6),".")
plot(time,vecnorm([traj(:,3),traj(:,6)]'),"--k")

figure
hold on
plot(traj(:,1),traj(:,4))
plot(xIMM(:,1),xIMM(:,3))

figure
hold on
plot([0 time],cBar(:,1))
plot([0 time],cBar(:,2))
plot([0 time],cBar(:,3))