clear;clc;close all

load("soloKF_CVCT.mat");
load("IMM_CVCT.mat");
load("trajectoryDataCVCT.mat");

figure 
subplot 211
plot(time(2:end),traj(2:end,1)-Ximm(1,:)',"linewidth",2)
hold on
grid on
plot(time(2:end),traj(2:end,1)-x1_hat(:,1),"linewidth",2)
plot(time(2:end),traj(2:end,1)-x2_hat(:,1),"linewidth",2)
plot(time(2:end),traj(2:end,1)-x3_hat(:,1),"linewidth",2)
xlabel("Time (sec)")
ylabel("X Error (m)")
legend(["IMM" "CV" "CT" "CA"],"location","best")
subplot 212
plot(time(2:end),traj(2:end,2)-Ximm(2,:)',"linewidth",2)
hold on 
grid on
plot(time(2:end),traj(2:end,2)-x1_hat(:,2),"linewidth",2)
plot(time(2:end),traj(2:end,2)-x2_hat(:,2),"linewidth",2)
plot(time(2:end),traj(2:end,2)-x3_hat(:,2),"linewidth",2)
xlabel("Time (sec)")
ylabel("Y Error (m)")

figure 
subplot 211
plot(time(2:end),traj(2:end,3)-Ximm(3,:)',"linewidth",2)
hold on
grid on
plot(time(2:end),traj(2:end,3)-x1_hat(:,3),"linewidth",2)
plot(time(2:end),traj(2:end,3)-x2_hat(:,3),"linewidth",2)
plot(time(2:end),traj(2:end,3)-x3_hat(:,3),"linewidth",2)
xlabel("Time (sec)")
ylabel("V_y Error (m)")
legend(["IMM" "CV" "CT" "CA"],"location","best")
subplot 212
plot(time(2:end),traj(2:end,4)-Ximm(4,:)',"linewidth",2)
hold on 
grid on
plot(time(2:end),traj(2:end,4)-x1_hat(:,4),"linewidth",2)
plot(time(2:end),traj(2:end,4)-x2_hat(:,4),"linewidth",2)
plot(time(2:end),traj(2:end,4)-x3_hat(:,4),"linewidth",2)
xlabel("Time (sec)")
ylabel("V_y Error (m)")

% figure 
% subplot 211
% plot(time(2:end),traj(2:end,1)-Ximm(1,:)',"linewidth",2)
% hold on
% grid on
% plot(time(2:end),traj(2:end,1)-x2_hat(:,1),"linewidth",2)
% xlabel("Time (sec)")
% ylabel("X Error (m)")
% legend(["IMM" "CT"],"location","best")
% subplot 212
% plot(time(2:end),traj(2:end,2)-Ximm(2,:)',"linewidth",2)
% hold on 
% grid on
% plot(time(2:end),traj(2:end,2)-x2_hat(:,2),"linewidth",2)
% xlabel("Time (sec)")
% ylabel("Y Error (m)")
% 
% figure 
% subplot 211
% plot(time(2:end),traj(2:end,3)-Ximm(3,:)',"linewidth",2)
% hold on
% grid on
% plot(time(2:end),traj(2:end,3)-x2_hat(:,3),"linewidth",2)
% xlabel("Time (sec)")
% ylabel("V_y Error (m)")
% legend(["IMM" "CT"],"location","best")
% subplot 212
% plot(time(2:end),traj(2:end,4)-Ximm(4,:)',"linewidth",2)
% hold on 
% grid on
% plot(time(2:end),traj(2:end,4)-x2_hat(:,4),"linewidth",2)
% xlabel("Time (sec)")
% ylabel("V_y Error (m)")