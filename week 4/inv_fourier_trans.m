    clc;
clear all;

% Read the spectral data from Excel files
heave_data = readmatrix('heave_amp_spec.xlsx');
pitch_data = readmatrix('pitch_amp_spec.xlsx');

% Extract wavenumber and amplitude spectrum
heave_wavenumber = heave_data(:,1);
heave_spectrum = heave_data(:,2);

pitch_wavenumber = pitch_data(:,1);
pitch_spectrum = pitch_data(:,2);

% Define time vector
num_samples = 100000; % Number of time samples
time = linspace(1000, 10000, num_samples); % Time range from 1000 to 10000

% Compute the inverse transform manually using summation
heave_time_series = zeros(size(time));
pitch_time_series = zeros(size(time));

for k = 1:length(heave_wavenumber)-1
    wave_angular_freq = sqrt(heave_wavenumber(k));
    heave_time_series = heave_time_series + heave_spectrum(k) * cos(wave_angular_freq * time);
end

for k = 1:length(pitch_wavenumber)-1
    wave_angular_freq = sqrt(pitch_wavenumber(k));
    pitch_time_series = pitch_time_series + pitch_spectrum(k) * cos(wave_angular_freq * time);
end

% Plot the results
figure;
subplot(2,1,1);
plot(time, heave_time_series, 'b');
xlabel('Time'); ylabel('Heave Response');
title('Heave Time Series'); grid on;

subplot(2,1,2);
plot(time, pitch_time_series, 'r');
xlabel('Time'); ylabel('Pitch Response');
title('Pitch Time Series'); grid on;

% Save the results to new Excel files
heave_output = [time' heave_time_series'];
pitch_output = [time' pitch_time_series'];

writematrix(heave_output, 'heave_time_series.xlsx');
writematrix(pitch_output, 'pitch_time_series.xlsx');

disp('Time series data has been saved to Excel files and plotted.');