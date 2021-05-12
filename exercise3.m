% In automotive communications, assume a bandwidth of BW = 500 MHz, a 
% carrier of fc = 80 GHz, a pulse duration of T = 40µs, and a range of
% detection of d = 300 m. Assume the target is moving at v = 150 km/h.
% Consider that we need to cover a range swath of 10m, so that you can 
% choose L and M = 4. Plot equation (9.34) to locate the Doppler shift
% and range.
%
% Eq. 9.34:
% l(delta_tau) <= 1/(pi^2*BW^2*delta_tau^2) 
%
% Luiz Felipe da S. Coelho  - luizfelipe.coelho@smt.ufrj.br
% may 2021
%

clc
clear
close all

% -------------------------------------------------------------------------
%                             Definitions
% -------------------------------------------------------------------------
kmh2ms = @(x) 1000*x/(60*60);  % Function to convert km/h to m/s
ms2kmh = @(x) (60*60)*x/1000;  % Function to convert m/s to km/h
c = 299792458;  % Speed of light, m/s
BW = 500*1e6;  % Bandwidth
f_c = 80*1e9;  % Carrier frequency
T = 40*1e-6;  % Pulse duration
dc = .2;  % Duty cycle
T_ch = dc*T;  % Chirp duration
d = 300;  % Range of detection
v = 150;  % Target's velocity
d_w = 10;  % Range swath
delta_d = c/(2*BW);  % Range resolution
L = round(d_w/delta_d);  % Total number of samples in fast-time
M = 4;  % Number of pulses
f_D = ((2*kmh2ms(v))/(c-kmh2ms(v)))*f_c;  % Doppler effect
nu = f_D/f_c;  % Doppler ratio
f_c2 = 0;  % Reduce computation burden.
a = BW/T;
Ts = kmh2ms(abs(v))*T*(1/c);
Fs = 1/Ts;
tau = 2*d/c;
n0 = round(tau*Fs - L/2);
fft_size = 2^18;
delta_d = linspace(n0, n0+L, L).*(c/(2*Fs));


% -------------------------------------------------------------------------
%                    Signal generation and processing
% -------------------------------------------------------------------------
% Time axis
t_ch = linspace(0, T_ch, T_ch*Fs);

% Transmitted signal (reference signal for xcorr)
tx = exp(1j*pi*(f_c2).*t_ch).*exp(1j*pi*a.*t_ch.^2);
tx = [tx zeros(1, round((T-T_ch)*Fs))];
ll = linspace(n0, n0+L, L);

% Received signal
target = exp(1j*pi*(f_c2+f_D).*t_ch).*exp(1j*pi*a*(1+nu).*t_ch.^2);
interf = exp(1j*pi*(f_c2+f_D).*t_ch).*exp(1j*pi*a*(1+nu).*t_ch.^2);

% Memory allocation
A = zeros(2*L + 1, M, L);
N = L+n0;
for m = 1:M
    
    d_m = d + (m-1)*kmh2ms(v)*T;
    t_samp_m = round(2*(d_m/c)*Fs);
    for d_idx = 1:L
        interf_samp = round(2*(delta_d(d_idx)/c)*Fs);
        target_m = [zeros(1, t_samp_m) target zeros(1, round((T-T_ch)*Fs) - t_samp_m)]*(2/3);
        interf_m = [zeros(1, interf_samp) interf zeros(1, round((T-T_ch)*Fs) - interf_samp)]*(1/3);
        rx_m = target_m + interf_m;
        % Dechirping process
        A(:, m, d_idx) = xcorr(rx_m(1, n0+1:N), tx(1, n0+1:N), L);
    end
end

A1 = A(L+2:end, :, 1);
LL = linspace(0, L-1, L);
MM = linspace(0, M-1, M);


figure(1),
surf(MM, LL, abs(A1))
xlabel('Slow-time, $m$', 'interpreter', 'latex')
ylabel('Fast-time, $l$', 'interpreter', 'latex')
zlabel('Normalized dechirped signal', 'interpreter', 'latex')
title('Slow-time vs. Fast-time matrix', 'interpreter', 'latex')
xticks([0 1 2 3])


A5 = A(L+2:end, :, 5);

figure(1),
surf(MM, LL, abs(A5))
xlabel('Slow-time, $m$', 'interpreter', 'latex')
ylabel('Fast-time, $l$', 'interpreter', 'latex')
zlabel('Normalized dechirped signal', 'interpreter', 'latex')
title('Slow-time vs. Fast-time matrix', 'interpreter', 'latex')
xticks([0 1 2 3])

AA = fftshift(fft(A1(:, 1), fft_size));
freq = fftshift([linspace(0, fft_size/2-1, fft_size/2)...
                 linspace(-fft_size/2, -1, fft_size/2)])*(Fs/fft_size);

figure(2),
plot(freq, abs(AA), 'linewidth', 2), hold on
% plot(freq(f_idx), val, 'x', 'markersize', 12, 'linewidth', 2), hold off, grid on
xlabel('Frequency, $f$ [Hz]', 'interpreter', 'latex')
ylabel('Magnitude', 'interpreter', 'latex')
title('Fast-time spectrum for $m=3$', 'interpreter', 'latex')
legend('Fast-time DFT', 'Maximum Value')
xlim([-2 2]*1e9)


delta_tau = linspace(0, 1, 200);
l = 1./(pi^2*BW^2.*delta_tau.^2);

figure(3),
plot(delta_tau, l, 'linewidth', 2), grid on
xlabel('$\Delta \tau$', 'interpreter', 'latex')
ylabel('$l(\Delta \tau)$', 'interpreter', 'latex')


% EoF

