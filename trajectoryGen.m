x0 = 0;
y0 = 0;
Vx0 = 1;
Vy0 = 0;
ax0 = 0;
ay0 = 0;


tf = 200;
fs   = 100;
dt   = 1/fs;

X = [x0;Vx0;ax0;y0;Vy0;ay0];

%% Straight Flight Mode

A = [1 dt .5*dt^2 0  0       0;
     0  1      dt 0  0       0;
     0  0       0 0  0       0;
     0  0       0 1 dt .5*dt^2;
     0  0       0 0  1      dt;
     0  0       0 0  0       0];

for ii = 1 : 100*fs

    X(:,ii+1) = A*X(:,ii);
    
end

%% Constant Angular Velocity Mode
w = 0.1;
A = [1      sin(w*dt)/w 0 0 -(1-cos(w*dt))/w 0;
     0        cos(w*dt) 0 0       -sin(w*dt) 0;
     0                0 0 0                0 0;
     0  (1-cos(w*dt))/w 0 1      sin(w*dt)/w 0;
     0        sin(w*dt) 0 0        cos(w*dt) 0;
     0                0 0 0                0 0]

for jj = ii : 150*fs

    X(:,jj+1) = A*X(:,jj);
    
end

%% Straight Flight Mode

A = [1 dt .5*dt^2 0  0       0;
     0  1      dt 0  0       0;
     0  0       1 0  0       0;
     0  0       0 1 dt .5*dt^2;
     0  0       0 0  1     dt;
     0  0       0 0  0       1];

for ii = jj : tf*fs

    X(:,ii+1) = A*X(:,ii) + [0 0 1 0 0 1]';
    
end

plot(X(1,:),X(4,:))
grid on
axis equal

save("trajectory","X")
