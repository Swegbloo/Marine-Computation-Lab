% MATLAB script for Pierson-Moskowitz spectrum, wave profile, animation, and FFT analysis
clc;
% Parameters
g = 9.81; % Acceleration due to gravity (m/s^2)
alpha = 8.1e-3; % Empirical constant
beta = 0.74; % Empirical constant
U = 10; % Wind speed (m/s)

% Frequency range (rad/s)
load('omega.mat'); % Avoid omega = 0 to prevent division by zero
domega = omega(2) - omega(1); % Frequency increment

% Pierson-Moskowitz spectrum formula
S_omega = (alpha * g^2 ./ omega.^5) .* exp(-beta * (g ./ (omega .* U)).^4);

% Compute amplitude spectrum
a_i = sqrt(2 * S_omega * domega);

% Compute wavenumber k using shallow water approximation
k = omega.^2 / g;

% Generate random phase angles (epsilon) for each frequency
epsilon = rand(size(omega)) * 2 * pi; % Random values between 0 and 2*pi

% Spatial and temporal grid for wave profile
x = linspace(0, 100, 500); % Spatial range (m)
t = linspace(0, 1000, 20000); % Temporal range (s) for higher frequency resolution

% Construct \eta(x,t) as a summation
[XX, TT] = meshgrid(x, t); % Create meshgrid for x and t
eta = zeros(size(XX)); % Initialize eta
for i = 1:length(omega)
    eta = eta + a_i(i) * cos(k(i) * XX - omega(i) * TT + epsilon(i));
end

% Task 1: Plot the Pierson-Moskowitz Spectrum and Amplitude Spectrum
figure;
subplot(2, 1, 1);
plot(omega, S_omega, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('\omega (rad/s)', 'FontSize', 12);
ylabel('S(\omega) (m^2s)', 'FontSize', 12);
title('Pierson-Moskowitz Spectrum: S(\omega) vs \omega', 'FontSize', 14);

subplot(2, 1, 2);
plot(omega, a_i, 'r-', 'LineWidth', 1.5);
grid on;
xlabel('\omega (rad/s)', 'FontSize', 12);
ylabel('a_i(\omega) (m)', 'FontSize', 12);
title('Amplitude Spectrum: a_i(\omega) vs \omega', 'FontSize', 14);