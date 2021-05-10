% Exercise 2:
%
% In automotive communications, the carrier frequency range lies between 76
% and 81 GHz. Assume a bandwidth of BW = 500 MHz, a carrier of f_c = 80
% GHz, a pulse duration of T = 40 \mus, and a range of detection d = 300 m.
% Assume the target is moving at v = 300 km/h. Consider that we need to
% cover a range swath of 20 m, so that you can choose L.
% 
% In this case, the Doppler effect originates a frequency shift of f_D =
% f_c*(2*v/c) = 44444.44 Hz.
% 
% Use a pulse burst FMCW with one and four pulses and estimate the range
% and the Doppler shift in both cases. Use a DFT in the fast-time axis to
% estimate the range, and then calculate the Doppler shift from the data
% delay for M = 1 and by exploring the slow-time information for M = 4.
%
% You can automatically detect a single target's delay by comparing the
% outcome of the matched filter with some prescribed threshold. Compare the
% results with those you know precisely the reflection delay, with the
% estimated delay, and with the exact delay +/- 5% of its nominal value.
%
%
% Luiz Felipe da S. Coelho - luizfelipe.coelho@smt.ufrj.br
% may 2021
%

clc
clear
close all

% -------------------------------------------------------------------------
%                               Definitions
% -------------------------------------------------------------------------
T = 40*1e-6;  % Pulse duration
dc = .2;  % Duty-cycle
T_ch = dc*T;  % Chirp duration
BW = 500*1e6;  % Signal bandwidth
f_c = 80*1e9;  % Carrier frequency
v = 300;  % Target speed, km/h
d = 300;  % Distance in meters
c = 299792458;  % Speed of light, m/s
d_w = 20;  % Range swath
delta_d = c/(2*BW);  % Range resolution
L = round(d_w/delta_d); % Total number of samples in fast-time
kmh2ms = @(x) 1000*x/(60*60);  % Function to convert km/h to m/s
f_D = ((2*kmh2ms(v))/(c-kmh2ms(v))) * f_c;  % Doppler effect
nu = f_D/f_c;  % Doppler ratio
f_c2 = 0;  % Reduce computational burden.
a = BW/T;  % Chirp rate
%
% Sampling rate estimation:
% To be able to detect a velocity of 300 km/h, one must have a minimal
% sampling rate. The sampling period is given by the following equation.
% Ts = v*(360/1000)*T*(1/c).
% To alleviate the high sampling rate, one can multiply the this time by M
% and use the 1st and the M-th pulse to estimate the velocity.
%
Ts = kmh2ms(abs(v))*T*(1/c);
Fs = 1/Ts;  % Sampling rate
M = 5;  % Number of pulses
d_est = 299.96;  % Arbitrary distance (what distance I wanna see)
n0 = round((2*d_est*Fs)/c);
fft_size = 2^18;


% -------------------------------------------------------------------------
%                    Signal generation and processing
% -------------------------------------------------------------------------
% Time axis
t_ch = linspace(0, T_ch, T_ch*Fs);

% Transmitted signal (reference signal for xcorr)
tx = exp(1j*pi*(f_c2).*t_ch).*exp(1j*pi*a.*t_ch.^2);
tx = [tx zeros(1, round((T-T_ch)*Fs))];

% Received signal
rx = exp(1j*pi*(f_c2+f_D).*t_ch).*exp(1j*pi*a*(1+nu).*t_ch.^2);

% Memory allocation
A = zeros(2*L + 1, M);
N = L+n0;
for m = 1:M
    d_m = d + (m-1)*kmh2ms(v)*T;
    t_samp_m = round(2*(d_m/c)*Fs);
    rx_m = [zeros(1, t_samp_m) rx zeros(1, round((T-T_ch)*Fs) - t_samp_m)];
    % Dechirping process
    A(:, m) = xcorr(rx_m(1, n0+1:N), tx(1, n0+1:N), L, 'normalized');
end

% -------------------------------------------------------------------------
%                   Slow-time vs. Fast-time matrix
% -------------------------------------------------------------------------
A = A(L+2:end, :);
LL = linspace(0, L-1, L);
MM = linspace(0, M-1, M);

figure,
surf(MM, LL, abs(A))
xlabel('Slow-time, $m$', 'interpreter', 'latex')
ylabel('Fast-time, $l$', 'interpreter', 'latex')
zlabel('Normalized dechirped signal', 'interpreter', 'latex')
title('Slow-time vs. Fast-time matrix', 'interpreter', 'latex')


% -------------------------------------------------------------------------
%                           Fast-time DFT
% -------------------------------------------------------------------------
fprintf('Distance estimation, (DFT):\n')
% M = 0
AA = fftshift(fft(A(:, 1), fft_size));
freq = fftshift([linspace(0, fft_size/2-1, fft_size/2)...
                 linspace(-fft_size/2, -1, fft_size/2)])*(Fs/fft_size);

[~, f_idx] = max(abs(AA));
f_t = freq(f_idx);
d_fft = (f_t*c)/(2*a);
fprintf('Estimeated range, for m = 0: %.4f m.\n', d_fft)

% M = 1
AA = fftshift(fft(A(:, 2), fft_size));
freq = fftshift([linspace(0, fft_size/2-1, fft_size/2)...
                 linspace(-fft_size/2, -1, fft_size/2)])*(Fs/fft_size);

[~, f_idx] = max(abs(AA));
f_t = freq(f_idx);
d_fft = (f_t*c)/(2*a);
fprintf('Estimeated range, for m = 1: %.4f m.\n', d_fft)

% M = 2
AA = fftshift(fft(A(:, 3), fft_size));
freq = fftshift([linspace(0, fft_size/2-1, fft_size/2)...
                 linspace(-fft_size/2, -1, fft_size/2)])*(Fs/fft_size);

[~, f_idx] = max(abs(AA));
f_t = freq(f_idx);
d_fft = (f_t*c)/(2*a);
fprintf('Estimeated range, for m = 2: %.4f m.\n', d_fft)

% M = 3
AA = fftshift(fft(A(:, 4), fft_size));
freq = fftshift([linspace(0, fft_size/2-1, fft_size/2)...
                 linspace(-fft_size/2, -1, fft_size/2)])*(Fs/fft_size);

[val, f_idx] = max(abs(AA));
f_t = freq(f_idx);
d_fft = (f_t*c)/(2*a);
fprintf('Estimeated range, for m = 3: %.4f m.\n', d_fft)
fprintf('\n')

figure,
plot(freq, abs(AA), 'linewidth', 2), hold on
plot(freq(f_idx), val, 'x', 'markersize', 12, 'linewidth', 2), hold off, grid on
xlabel('Frequency, $f$ [Hz]', 'interpreter', 'latex')
ylabel('Magnitude', 'interpreter', 'latex')
title('Fast-time spectrum for $m=3$', 'interpreter', 'latex')
legend('Fast-time DFT', 'Maximum Value')
xlim([-2 2]*1e9)

% -------------------------------------------------------------------------
%                           Distance Estiamation
% -------------------------------------------------------------------------
fprintf('Distance estimation (maximum value in matrix):\n')
% m = 0
[~, idx] = max(abs(A(:, 1)));
d_hat1 = ((LL(idx)+n0)/Fs)*c/2;
fprintf('Estimated distance m = 0: %.4f m.\n', d_hat1)
% m = 1
[~, idx] = max(abs(A(:, 2)));
d_hat2 = ((LL(idx)+n0)/Fs)*c/2;
fprintf('Estimated distance m = 1: %.4f m.\n', d_hat2)
% m = 2
[~, idx] = max(abs(A(:, 3)));
d_hat3 = ((LL(idx)+n0)/Fs)*c/2;
fprintf('Estimated distance m = 2: %.4f m.\n', d_hat3)
% m = 3
[val, idx] = max(abs(A(:, 4)));
d_hat4 = ((LL(idx)+n0)/Fs)*c/2;
fprintf('Estimated distance m = 3: %.4f m.\n', d_hat4)
fprintf('\n')


% -------------------------------------------------------------------------
%                           Velocity Estimation
% -------------------------------------------------------------------------
fprintf('Velocity estimation:\n')
vel = (d_hat2-d_hat1)/(T);
ms2kmh = @(x) (60*60)*x/1000;  % Function to convert m/s to km/h
fprintf('Estimated velocity, using m = 0 and m = 1: %.2f km/h.\n',...
        ms2kmh(vel))

vel = (d_hat3-d_hat1)/(2*T);
ms2kmh = @(x) (60*60)*x/1000;  % Function to convert m/s to km/h
fprintf('Estimated velocity, using M = 0 and M = 2: %.2f km/h.\n',...
        ms2kmh(vel))

vel = (d_hat4-d_hat1)/(3*T);
ms2kmh = @(x) (60*60)*x/1000;  % Function to convert m/s to km/h
fprintf('Estimated velocity, using M = 0 and M = 3: %.2f km/h.\n',...
        ms2kmh(vel))
    
% EoF
