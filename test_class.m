% MATLAB script for Pierson-Moskowitz spectrum, wave profile, animation, and FFT analysis
clc;
% Parameters
g = 9.81; % Acceleration due to gravity (m/s^2)
alpha = 8.1e-3; % Empirical constant
beta = 0.74; % Empirical constant
U = 10; % Wind speed (m/s)

% Frequency range (rad/s)
omega = linspace(0.1, 3, 1000); % Avoid omega = 0 to prevent division by zero
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

% Task 2 & 3: Movie of irregular wave profile
figure;
for n = 1:200:length(t)
    plot(x, eta(n, :), 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('x (m)', 'FontSize', 12);
    ylabel('\eta(x,t) (m)', 'FontSize', 12);
    title(sprintf('Irregular Wave Profile: \\eta(x,t) at t = %.2f s', t(n)), 'FontSize', 14);
    pause(0.05); % Adjust for smoother animation
end

% Task 4: FFT Analysis with Minimization
% Choose a fixed x and generate time history
x_fixed = 50; % Fixed position along x
eta_time_history = interp2(XX, TT, eta, x_fixed, t); % Extract time history

% Plot the time history of eta at x_fixed
figure;
plot(t, eta_time_history, 'k-', 'LineWidth', 1);
grid on;
xlabel('Time (s)', 'FontSize', 12);
ylabel('\eta(x=50, t) (m)', 'FontSize', 12);
title('Time History of Surface Elevation at x = 50 m', 'FontSize', 14);

% Perform FFT on the time history
N = length(t);
eta_fft = fft(eta_time_history) / N; % Normalize FFT by the number of points
freq = (0:N/2-1) * (1 / (t(end) - t(1))); % Frequency axis for one-sided FFT
eta_fft_magnitude = abs(eta_fft(1:N/2)); % One-sided FFT magnitude

% Map FFT frequency axis to the original omega grid
fft_interp = interp1(2 * pi * freq, eta_fft_magnitude, omega, 'linear', 'extrap');

% Apply higher-order smoothing (5-point moving average)
fft_smoothed = movmean(fft_interp, 5);

% Optimize scaling factor to minimize error
gain_factor = sum(a_i .* fft_smoothed) / sum(fft_smoothed.^2);
fft_optimized = gain_factor * fft_smoothed;

% Compute the area under the FFT signal and the original amplitude spectrum
area_fft = trapz(omega, fft_optimized);
area_original = trapz(omega, a_i);

% Compare Optimized FFT result with amplitude spectrum
figure;
plot(omega, fft_optimized, 'b-', 'LineWidth', 1.5); % Optimized FFT magnitude
hold on;
plot(omega, a_i, 'r--', 'LineWidth', 1.5); % Original amplitude spectrum
grid on;
xlabel('\omega (rad/s)', 'FontSize', 12);
ylabel('Amplitude (m)', 'FontSize', 12);
title(sprintf('Optimized FFT vs Original Amplitude Spectrum\nArea FFT: %.4f, Area Original: %.4f', area_fft, area_original), 'FontSize', 14);
legend('Optimized FFT of Time History', 'Original Amplitude Spectrum');