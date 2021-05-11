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
d = 300;  % Range of detection
v = 150;  % Target's velocity
d_w = 10;  % Range swath
delta_d = c/(2*BW);  % Range resolution
L = round(d_w/delta_d);  % Total number of samples in fast-time
M = 4;  % Number of pulses

delta_tau = linspace(0, 1, 200);
l = 1./(pi^2*BW^2.*delta_tau.^2);

figure(1),
plot(delta_tau, l, 'linewidth', 2), grid on
xlabel('$\Delta \tau$', 'interpreter', 'latex')
ylabel('$l(\Delta \tau)$', 'interpreter', 'latex')