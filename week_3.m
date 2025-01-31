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

% Extract time history at x_fixed
x_fixed = 50; % Fixed position along x
eta_time_history = interp2(XX, TT, eta, x_fixed, t); % Extract time history

% Significant Wave Height Calculation
% Zero Up-Crossing Method
zero_crossings_up = find(diff(sign(eta_time_history)) > 0);
wave_heights_up = abs(diff(eta_time_history(zero_crossings_up)));
Hs_up = 4 * std(wave_heights_up); % Hs estimation

% Zero Down-Crossing Method
zero_crossings_down = find(diff(sign(eta_time_history)) < 0);
wave_heights_down = abs(diff(eta_time_history(zero_crossings_down)));
Hs_down = 4 * std(wave_heights_down); % Hs estimation

% Peak Wave Period (Tp)
eta_fft = fft(eta_time_history) / length(t);
freq = (0:length(t)/2-1) * (1 / (t(end) - t(1)));
[~, peak_idx] = max(abs(eta_fft(1:length(freq))));
Tp = 1 / freq(peak_idx); % Peak period

% Display results
fprintf('Significant Wave Height (Zero Up-Crossing): %.3f m\n', Hs_up);
fprintf('Significant Wave Height (Zero Down-Crossing): %.3f m\n', Hs_down);
fprintf('Peak Wave Period: %.3f s\n', Tp);

% Plot the time history of eta at x_fixed
figure;
plot(t, eta_time_history, 'k-', 'LineWidth', 1);
grid on;
xlabel('Time (s)', 'FontSize', 12);
ylabel('\eta(x=50, t) (m)', 'FontSize', 12);
title('Time History of Surface Elevation at x = 50 m', 'FontSize', 14);
