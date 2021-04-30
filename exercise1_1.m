% • For a chirp supported for a normalized period of time Tch = 1 = T plot
% the real part of the radar chirp for aT = 20.
% • Plot the ambiguity function.
% • Plot the spectrum of the radar signal.
% • Repeat the above plots for sinusoid pulse and for a square wave pulse.

clc
clear
close all

Tch = 1;
T = 1;
a = 20;
Fs = 1e3;
t = linspace(0, T, T*Fs);
f_c = 0;

chirp = exp(1j*pi*(f_c).*t).*exp(1j*(pi*a.*t.^2));
[r_chirp, lags] = xcorr(chirp, 10, 'normalized');
fft_size = 2^ceil(log2(length(chirp)));

CHIRP = fftshift(fft(chirp, fft_size))/(fft_size);
freq = linspace(-0.5, 0.5, fft_size)/Fs;

figure(1),
plot(t, real(chirp)), grid on
xl=xlabel('Time, $t$ [sec]', 'interpreter', 'latex')
%ylabel('Amplitude, $x_{\textrm{chirp}}(t)$', 'interpreter', 'latex')
yl=ylabel('Amplitude, $x(t)$', 'interpreter', 'latex')
title('Real part of complex chirp signal', 'interpreter', 'latex')
set(xl,'Interpreter','latex');set(xl,'FontSize',14);
set(yl,'FontSize',14);

zoom on;
figure(2),
plot(lags, abs(r_chirp)), grid on
xl=xlabel('Lags [samples]', 'interpreter', 'latex')
%ylabel('$R_{xx}$', 'interpreter', 'latex')
yl=ylabel('$A(t,f_{\rm D})$', 'interpreter', 'latex')
title('Ambiguity Function', 'interpreter', 'latex')
set(xl,'Interpreter','latex');set(xl,'FontSize',14);
set(yl,'FontSize',14);

figure(3),
plot(freq,abs(CHIRP)), grid on
xl=xlabel('Frequency, $f$ [Hz]', 'interpreter', 'latex');
yl=ylabel('Magnitude $|X(f)|$', 'interpreter', 'latex');
title('Spectrum of the Radar Signal', 'interpreter', 'latex')
set(xl,'Interpreter','latex');set(xl,'FontSize',14);
set(yl,'FontSize',14);
xlim([-1.5*10^(-4) 1.5*10^(-4)]);


% Sinusoid pulse
f_c = 10;  % sinusoid frequency
phi = 0;  % phase

sinusoid = exp(1j*(2*pi*f_c*t + phi));
[r_sin, lags] = xcorr(sinusoid, 10, 'normalized');
fft_size = 2^ceil(log2(length(sinusoid)));

SINUSOID = fftshift(fft(sinusoid, fft_size))/(fft_size);
freq = linspace(-.5, .5, fft_size)/Fs;

figure(4),
plot(t, real(sinusoid)), grid on
xl=xlabel('Time, $t$ [sec]', 'interpreter', 'latex')
yl=ylabel('Amplitude, $x(t)$', 'interpreter', 'latex')
%ylabel('Amplitude, $x_{\textrm{chirp}}(t)$', 'interpreter', 'latex')
title('Real part of complex sinusoidal signal', 'interpreter', 'latex')
set(xl,'Interpreter','latex');set(xl,'FontSize',14);
set(yl,'FontSize',14);
%ylim([-0.5 1.5]);

figure(5),
plot(lags, abs(r_sin)), grid on
xl=xlabel('Lags [samples]', 'interpreter', 'latex')
%ylabel('$R_{xx}$', 'interpreter', 'latex')
yl=ylabel('$A(t,f_{\rm D})$', 'interpreter', 'latex')
title('Ambiguity Function', 'interpreter', 'latex')
set(xl,'Interpreter','latex');set(xl,'FontSize',14);
set(yl,'FontSize',14);


figure(6),
plot(freq, abs(SINUSOID)), grid on
xl=xlabel('Frequency, $f$ [Hz]', 'interpreter', 'latex')
yl=ylabel('Magnitude $|X(f)|$', 'interpreter', 'latex')
title('Spectrum of the Radar Signal', 'interpreter', 'latex')
set(xl,'Interpreter','latex');set(xl,'FontSize',14);
set(yl,'FontSize',14);
xlim([-1.5*10^(-4) 1.5*10^(-4)]);

% Square wave
d_c = .4;  % Duty-cycle
T = 1;  % Time duration
Fs = 1e3;  % Frequency sample

t = linspace(0, T, T*Fs);  % Time axis
sqr = [ones(1, (d_c*T*Fs)) zeros(1, ((1-d_c)*T*Fs))];
[r_sqr, lags] = xcorr(sqr, 10, 'normalized');
fft_size = 2^ceil(log2(length(sqr)));

SQR = fftshift(fft(sqr, fft_size))/(fft_size);
freq = linspace(-.1, .1, fft_size)/Fs;

figure(7),
plot(t, real(sqr)), grid on
xl=xlabel('Time, $t$ [sec]', 'interpreter', 'latex')
yl=ylabel('Amplitude, $x(t)$', 'interpreter', 'latex')
%ylabel('Amplitude, $x_{\textrm{square}}(t)$', 'interpreter', 'latex')
title('Real part of complex sinusoidal signal', 'interpreter', 'latex')
yl=ylabel('$A(t,f_{\rm D})$', 'interpreter', 'latex')
set(xl,'Interpreter','latex');set(xl,'FontSize',14);
set(yl,'FontSize',14);
ylim([-0.5 1.5]);

figure(8),
plot(lags, abs(r_sqr)), grid on
xl=xlabel('Lags [samples]', 'interpreter', 'latex')
%ylabel('$R_{xx}$', 'interpreter', 'latex')
yl=ylabel('$A(t,f_{\rm D})$', 'interpreter', 'latex')
title('Ambiguity Function', 'interpreter', 'latex')
set(xl,'Interpreter','latex');set(xl,'FontSize',14);
set(yl,'FontSize',14);


figure(9),
plot(freq, abs(SQR)), grid on
xl=xlabel('Frequency, $f$ [Hz]', 'interpreter', 'latex')
yl=ylabel('Magnitude $|X(f)|$', 'interpreter', 'latex')
title('Spectrum of the Radar Signal', 'interpreter', 'latex')
set(xl,'Interpreter','latex');set(xl,'FontSize',14);
set(yl,'FontSize',14);
xlim([-1*10^(-5) 1*10^(-5)]);



