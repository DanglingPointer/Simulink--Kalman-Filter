clear
syms K
syms T

%Vi har lagt til et sinussignal inn i systemet og fjernet de andre
%påvirkningene. Da kan vi se hva den rene påvikningne fra bølgene er i pent
%vær. 

w_1 = 0.005;    %Frekvensen til Sin_1
w_2 = 0.05;     %Frekvensen til Sin_2
a_1 = 31.9755; % 0.5581;   %Amplituden til Sin_1
a_2 = 0.7775; %0.0136;   %Amplituden til Sin_2

eqv_1 = (K/T)/sqrt(w_1^4 +w_1^2/T^2);  %Transferfunksjonen til Sin_1

eqv_2 = (K/T)/sqrt(w_2^4 +w_2^2/T^2);  %Transferfunksjonen til Sin_2

[K, T] = vpasolve([eqv_1 == a_1, eqv_2 == a_2], [K, T]); %Lager vektor og regner ut K og T verdi

K = double(K); % K = 0.1745
T = double(T); % T = 87.5282

% Designing PD-controller
omega_c = 0.10; % [rad/s]
T_f = -1/(omega_c*tan(130*pi/180));         % T_f = 8.3910
K_pd = sqrt(T_f^2*omega_c^4+omega_c^2)/K;   % K_pd = 0.7480