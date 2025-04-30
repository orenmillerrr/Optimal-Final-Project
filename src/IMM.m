clear;clc;close all

load("trajectory_data.mat")

n = length(traj);
xySigma = 1;
VxySigma = .5;
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
H = [0 1 0 0;
     0 0 0 1];

x0 = [traj(1,1) traj(1,3) traj(1,2) traj(1,4)]';
P0 = eye(4);

kf1 = kalmanFilter(A,[],H,Q,R,x0,P0);

%% Kalman Filter 2
qTURN = 25;
Q = qCV*dt*[1 0;0 1];
Q = [Q zeros(2);zeros(2) Q];

w = deg2rad(90/5);
A = [1    sin(w*dt)/w    0  -(1-cos(w*dt))/w;
      0      cos(w*dt)   0        -sin(w*dt);
      0 (1-cos(w*dt))/w  1      sin(w*dt)/w;
      0       sin(w*dt)  0         cos(w*dt)];

kf2 = kalmanFilter(A,[],H,Q,R,x0,P0);

%% IMM Parameters
x{1}  = x0;
P{1}  = P0;
x{2}  = x0;
P{2}  = P0;
M = [0.97 0.03;0.03 0.97];    % Markov Chain
mu = [0.5 0.5];
cBar = mu * M;
numModels = 2;

for ii = 1 : length(time)

    [Z1,S1] = kf1.computeInnovation(y(ii,:)');
    [L1] = kf1.computeLikelihood(Z1,S1);

    [Z2,S2] = kf2.computeInnovation(y(ii,:)');
    [L2] = kf2.computeLikelihood(Z2,S2);

    % IMM Likelyhood
    likeIMM = [L1 L2];

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

    xIMM(ii,:) = weightsIMM(1)*x{1} + weightsIMM(2)*x{2};
    pIMM = (P{1} + (x{1} -  xIMM(ii,:))'*(x{1} -  xIMM(ii,:)))*weightsIMM(1) + ...
           (P{2} + (x{2} -  xIMM(ii,:))'*(x{2} -  xIMM(ii,:)))*weightsIMM(2);

    % Compute Mixed Initial Conditions
    cBar1= M(:,1)'*weightsIMM';  
    if cBar1>1e-80
        omega(1,1)=M(1,1)*weightsIMM(1)/cBar1;
        omega(2,1)=M(2,1)*weightsIMM(2)/cBar1;
    else
        cBar1 = 0;
        omega(1,1) = 0;
        omega(2,1) = 0;
    end
    cBar2=M(:,2)'*weightsIMM';
    if cBar2>1e-80
        omega(1,2)=M(1,2)*weightsIMM(1)/cBar2;
        omega(2,2)=M(2,2)*weightsIMM(2)/cBar2;
    else
        cBar2 = 0;
        omega(1,2) = 0;
        omega(2,2) = 0;
    end
    cBar(ii+1,:)=[cBar1,cBar2]; 

    x1Mixed = omega(1,1)*x{1}+omega(2,1)*x{2};
    x2Mixed = omega(1,2)*x{1}+omega(2,2)*x{2};

    P1Mixed = (kf1.P + (kf1.x_hat - x1Mixed)'*(kf1.x_hat - x1Mixed))*omega(1,1) + ...
              (kf2.P + (kf2.x_hat - x2Mixed)'*(kf2.x_hat - x2Mixed))*omega(2,1);

    P2Mixed = (kf1.P + (kf1.x_hat - x1Mixed)'*(kf1.x_hat - x1Mixed))*omega(1,2) + ...
              (kf2.P + (kf2.x_hat - x2Mixed)'*(kf2.x_hat - x2Mixed))*omega(2,2);

    kf1.x_hat = x1Mixed;
    kf2.P = P1Mixed;

    kf2.x_hat = x2Mixed;
    kf2.P = P1Mixed;

    % Time Update
    kf1.propagate
    kf2.propagate

end

figure
hold on
plot(traj(:,1),traj(:,2))
plot(xIMM(:,1),xIMM(:,3))

figure
plot(cBar)