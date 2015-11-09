function [ psd ] = PowerDensity( lambda, omega )
[n_rows ,n_col] = size(omega);
psd = nan(n_rows, 1);

sigma = 0.1210;
omega_0 = 0.7823;

for i = 1:n_rows
    psd(i) = ((2*lambda*omega_0*sigma)^2*omega(i)^2)/...
    (omega(i)^4-2*omega_0^2*omega(i)^2+4*lambda^2*omega_0^2*omega(i)^2+omega_0^4);
end

end

