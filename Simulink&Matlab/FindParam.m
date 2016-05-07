function [ K,T, omega_0, lambda, K_w] = FindParam()
% ----------------------Kopiert fra oppg5_1_b------------------------------
syms K
syms T
w_1 = 0.005;    %Frekvensen til Sin_1
w_2 = 0.05;     %Frekvensen til Sin_2
a_1 = 31.9755; % 0.5581;   %Amplituden til Sin_1
a_2 = 0.7775; %0.0136;   %Amplituden til Sin_2
eqv_1 = (K/T)/sqrt(w_1^4 +w_1^2/T^2);  %Transferfunksjonen til Sin_1
eqv_2 = (K/T)/sqrt(w_2^4 +w_2^2/T^2);  %Transferfunksjonen til Sin_2
[K, T] = vpasolve([eqv_1 == a_1, eqv_2 == a_2], [K, T]); %Lager vektor og regner ut K og T verdi
K = double(K); % K = 0.1745
T = double(T); % T = 87.5282
% ----------------------Kopiert fra oppg5_2_d------------------------------
load('wave.mat');
x = psi_w(2,:);     % waves' influence on the measurements, i.e. psi_w [degree]
x = x * (pi/180);   % [degree] to [rad]
fs = 10;            % sampling frequency [Hz]
window = 4096;      % window size
noverlap = [];      % not to be specified
nfft = [];          % not to be specified
[psd, f] = pwelch(x,window,noverlap,nfft,fs);
psd = psd*(1/(2*pi));     % [power/Hz] to [power*s/rad]
omega = f*2*pi;         % [Hz] to [rad/s]
[value, index] = max(psd);
omega_0 = omega(index); % omega_0 = 0.7823
sigma = sqrt(value);
OPTIONS=optimset('lsqcurvefit');
OPTIONS.TolFun = 1e-18;
OPTIONS.TolX=1e-40;
OPTIONS.MaxFunEvals=1e5;
lambda = lsqcurvefit('PowerDensity',.1,omega,psd); % lambda = 0.0074
K_w = 2*lambda*sigma; % K_w = 5.7006e-04

end

