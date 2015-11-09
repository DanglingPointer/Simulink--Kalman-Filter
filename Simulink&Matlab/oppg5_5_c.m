clear
[ K,T, omega_0, lambda, K_w] = FindParam(); % obtaining values
A = [0 1 0 0 0; -(omega_0^2) -2*lambda 0 0 0; ...
    0 0 0 1 0; 0 0 0 -1/T -K/T; 0 0 0 0 0];

B = zeros(5,4);
B(:,1) = [0 0 0 K/T 0]';
B(2,2)=K_w;
B(5,3)=1;
C = [0 1 1 0 0];
D = [0 0 0 1];

% Using Matlab functions
Ts = 0.1; % sampling time
sys = ss(A,B,C,D); % continuous state-space model
sysd = c2d(sys, Ts); % using Zero-order hold by default, see 
                     % http://se.mathworks.com/help/control/ref/c2d.html?searchHighlight=c2d#inputarg_method
% Parameters for Kalman filter (see assignment 5 problem 2)
Ad = sysd.a;
Bd = sysd.b(:,1);
Gd = sysd.b(:, 2:3);