x0 = 0;
y0 = 0;
Vx0 = 100;
Vy0 = 0;


tf = 200;
fs   = 100;
dt   = 1/fs;

X = [x0;Vx0;y0;Vy0];

%% Straight Flight Mode

A = [1 dt 0 0;
     0  1 0 0;
     0  0 1 dt;
     0  0 0  1];

for ii = 1 : 100*fs

    X(:,ii+1) = A*X(:,ii);
    
end

%% Constant Angular Velocity Mode
w = 0.05;
A = [1    sin(w*dt)/w   0  -(1-cos(w*dt))/w;
     0      cos(w*dt)   0        -sin(w*dt);
     0 (1-cos(w*dt))/w   1      sin(w*dt)/w;
     0       sin(w*dt)  0         cos(w*dt)];

for jj = ii : 150*fs

    X(:,jj+1) = A*X(:,jj);
    
end

%% Straight Flight Mode

A = [1 dt 0 0;
     0  1 0 0;
     0  0 1 dt;
     0  0 0  1];

for ii = jj : tf*fs

    X(:,ii+1) = A*X(:,ii);
    
end

plot(X(1,:),X(3,:))
grid on
axis equal

save("trajectory","X")
