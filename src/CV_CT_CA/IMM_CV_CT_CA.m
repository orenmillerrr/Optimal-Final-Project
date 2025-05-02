clear;clc;close all

load("trajectory/trajectoryDataCVCTCA.mat")

n = length(traj);
xySigma  = .5;
VxySigma = .1;
AxySigma = .01;
noise = [xySigma*randn([2 n]);
         VxySigma*randn([2 n]);
         AxySigma*randn([2 n])];
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
        y = meas(:,1:4);
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
qa = 10;
Qcv = [qp*dt*eye(2)     zeros(2)     zeros(2);
           zeros(2) qv*dt*eye(2)     zeros(2);
           zeros(2)     zeros(2) qa*dt*eye(2);];

kfCV = kalmanFilter(Acv,[],H,Qcv,R,x0,P0);

%% Constant Turn Kalman Filter
w = deg2rad(30);
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


%% Constant Velocity Kalman Filter
Aca = [eye(2) dt*eye(2) (dt^2)/2*eye(2);
       zeros(2)    eye(2) dt*eye(2);
       zeros(2)  zeros(2) eye(2)];
qp = 1;
qv = 1;
qa = 1;
Qca = [qp*dt*eye(2)     zeros(2)     zeros(2);
           zeros(2) qv*dt*eye(2)     zeros(2);
           zeros(2)     zeros(2) qa*dt*eye(2);];

kfCa = kalmanFilter(Aca,[],H,Qca,R,x0,P0);

%% IMM Initialization
x{1}  = x0;
P{1}  = P0;
x{2}  = x0;
P{2}  = P0;
x{3}  = x0;
P{3}  = P0;
M = [0.980 0.010 0.010;
     0.010 0.980 0.010;
     0.010 0.010 0.980];    % Markov Chain
mu = [1/3 1/3 1/3];
numModels = 3;

KF = [kfCV,kfCT,kfCa];

for ii = 1 : length(time)-1
    
    for r = 1 : numModels
        Xj{r} = KF(r).x_hat;
        Pj{r} = KF(r).P;
    end

    %% State Interaction
    % Compute Mode Probability
    cBar = mu(ii,:)*M;

    % Compute Conditional Model Probabilities
    for i = 1 : numModels
        for j = 1 : numModels
            w(i,j) = (M(i,j) * mu(ii,i)) / cBar(j);
        end
    end

    % Compute Mixed State Estimate for Each Model
    X0j{1} = w(1,1)*Xj{1} + w(2,1)*Xj{2} + w(3,1)*Xj{3};
    X0j{2} = w(1,2)*Xj{1} + w(2,2)*Xj{2} + w(3,2)*Xj{3};
    X0j{3} = w(1,3)*Xj{1} + w(2,3)*Xj{2} + w(3,3)*Xj{3};

    P0j{1} = (Pj{1} + (Xj{1} - X0j{1})'*(Xj{1} - X0j{1}))*w(1,1) + ...
              (Pj{2} + (Xj{2} - X0j{2})'*(Xj{2} - X0j{2}))*w(2,1) + ...
              (Pj{3} + (Xj{3} - X0j{3})'*(Xj{3} - X0j{3}))*w(3,1);

    P0j{2} = (Pj{1} + (Xj{1} - X0j{1})'*(Xj{1} - X0j{1}))*w(1,2) + ...
              (Pj{2} + (Xj{2} - X0j{2})'*(Xj{2} - X0j{2}))*w(2,2) + ...
              (Pj{3} + (Xj{3} - X0j{3})'*(Xj{3} - X0j{3}))*w(3,3);

    P0j{3} = (Pj{1} + (Xj{1} - X0j{1})'*(Xj{1} - X0j{1}))*w(1,3) + ...
              (Pj{2} + (Xj{2} - X0j{2})'*(Xj{2} - X0j{2}))*w(2,3) + ...
              (Pj{3} + (Xj{3} - X0j{3})'*(Xj{3} - X0j{3}))*w(3,3);

    %% Update Kalman Filters
    kfCV.x_hat = X0j{1};
    kfCV.propagate
    kfCV.update(y(ii,:)')
    kfCT.x_hat = X0j{2};
    kfCT.propagate
    kfCT.update(y(ii,:)')
    kfCa.x_hat = X0j{3};
    kfCa.propagate
    kfCa.update(y(ii,:)')

    like = [kfCV.like,kfCT.like, kfCa.like];
    %% Model Probability Update
    mu(ii+1,:) = like.*cBar / sum(like.*cBar);

    for r = 1 : numModels
        Xj{r} = KF(r).x_hat;
        Pj{r} = KF(r).P;
    end

    %% State Estimation Combination
    Ximm(:,ii) = mu(ii+1,1)*Xj{1} + ...
                  mu(ii+1,2)*Xj{2} + ...
                  mu(ii+1,3)*Xj{3};

    Pimm{ii} = (Pj{1} + (Xj{1}-Ximm(:,ii))'*(Xj{1}-Ximm(:,ii)))*mu(ii+1,1) + ...
                (Pj{2} + (Xj{2}-Ximm(:,ii))'*(Xj{2}-Ximm(:,ii)))*mu(ii+1,2) + ...
                (Pj{3} + (Xj{3}-Ximm(:,ii))'*(Xj{3}-Ximm(:,ii)))*mu(ii+1,3);

end

fileName = 'D:\Classes\MECH 7710 - Optimal Control\Optimal-Final-Project\src\savedData\IMM_CVCTCA';
save(fileName,'Ximm','time')
red = [1,0.7,.7,1];
green = [.63,.87,.7,1];
yellow = [1,1,.7,1];
blue = [0.7,0.8,1,.5];

% Plot the trajectory
figure;
plot(traj(1:10*fs,1),traj(1:10*fs,2), 'Color','r','LineWidth',2);
hold on 
grid on;
plot(traj(10*fs:20*fs,1),traj(10*fs:20*fs,2), 'Color','b','LineWidth',2);
plot(traj(20*fs:30*fs,1),traj(20*fs:30*fs,2), 'Color','g','LineWidth',2);
xlabel('X Position');
ylabel('Y Position');
legend(["CV" "CT"],"Location","best")
axis equal;

% figure 
% subplot 211
% plot(time(2:end),traj(2:end,1)-Ximm(1,:)')
% subplot 212
% plot(time(2:end),traj(2:end,2)-Ximm(2,:)')


figure 
plot(time,mu(:,1),"LineWidth",2,"color","r")
grid on
hold on
plot(time,mu(:,2),"LineWidth",2,"color","b")
plot(time,mu(:,3),"LineWidth",2,"color","g")
legend(["CV" "CT" "CA"],"location","best")

