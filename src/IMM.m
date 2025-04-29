clear;clc;close all

traj = load("trajectory.mat");
traj = traj.X;

n = length(traj);
xySigma = 5;
VxySigma = .5;
noise = [xySigma*randn([1 n]);
         xySigma*randn([1 n]);
         VxySigma*randn([1 n]);
         VxySigma*randn([1 n])];
y = traj + noise;

grid on
hold on
axis equal

tf = 200;
fs = 100;
dt = 1/fs;
time = 0 : dt : tf;
R = diag([xySigma xySigma VxySigma VxySigma]);

%% Kalman Filter 1
qCV = 400;
Q = qCV*dt*[dt^2/3 dt/2;dt/2 1];
Q = [Q zeros(2);zeros(2) Q];
A = [1 dt 0 0;
      0  1 0 0;
      0  0 1 dt;
      0  0 0  1];
H = [0 0 0 0;
      0 1 0 0;
      0 0 0 0;
      0 0 0 1];

x0 = traj(:,1);
P0 = eye(4);

kf1 = kalmanFilter(A,[],H,Q,R,x0,P0);

%% Kalman Filter 2
qTURN = 25;
Q = qCV*dt*[1 0;0 1];
Q = [Q zeros(2);zeros(2) Q];

w = 0.05;
A = [1    sin(w*dt)/w   0  -(1-cos(w*dt))/w;
      0      cos(w*dt)   0        -sin(w*dt);
      0 (1-cos(w*dt))/w   1      sin(w*dt)/w;
      0       sin(w*dt)  0         cos(w*dt)];

kf2 = kalmanFilter(A,[],H,Q,R,x0,P0);

%% IMM Parameters
M = [0.97 0.03;0.05 0.95];    % Markov Chai
x{1}  = x0;
P{1}  = P0;
x{2}  = x0;
P{2}  = P0;
mu = [0.7 0.3];
cBar = mu * M;
numModels = 2;

for ii = 1 : length(time)
    
    %% Run Kalman Filter 1
    % Time Update
    kf1.propagate

    % Measurement Update
    kf1.update(y(:,ii))

    x1_hat(:,ii) = kf1.getState;
    x{1} = x1_hat(:,ii);
    P{1} = kf1.getCovariance;
    L{1} = kf1.computeLikelihood;

    
    %% Run Kalman Filter 2
    % Time Update
    kf2.propagate

    % Measurement Update
    kf2.update(y(:,ii))

    x2_hat(:,ii) = kf2.getState;
    x{2} = x2_hat(:,ii);
    P{2} = kf2.getCovariance;
    L{2} = kf2.computeLikelihood;

    X(:,ii) = (mu * [x{1},x{2}]')';

    mu = [norm(L{1}*cBar(1)) norm(L{2}*cBar(2))] / sum([norm(L{1}*cBar(1)) norm(L{2}*cBar(2))]);
    cBar = mu * M;
    omega = zeros(numModels,numModels);
    for i = 1 : numModels
        for j = 1 : numModels
            omega(i,j) = (M(i,j) * mu(i)) / cBar(j);
        end
    end

    for j = 1 : numModels
        xMixed{j} = zeros(size(x{j}));
        PMixed{j} = zeros(size(P{j}));
        for i = 1 : numModels
            xMixed{j} = xMixed{j} + (omega(i,j) * x{j});
            PMixed{j} = PMixed{j} + (omega(i,j) * (P{j} + (x{j} - xMixed{j})*(x{j} - xMixed{j})'));
        end
    end

    kf1.immUpdate(x{1},P{1})
    kf2.immUpdate(x{2},P{2})

end

figure
hold on
plot(traj(1,:),traj(3,:))
plot(x1_hat(1,:),x1_hat(3,:),".")
plot(x2_hat(1,:),x2_hat(3,:),".")
plot(X(1,:),X(3,:))