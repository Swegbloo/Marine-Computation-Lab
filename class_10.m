clc
d =10; %depth
T = 3; %Time period
g = 9.81; % gravity
w = (2*pi)/T; %frequency
A = 0.1; %amplitude

disp(w);

k = 0.4473;
x = linspace(0, 20, 500); % x-axis range
t_interval = 0:0.1:10;     % Time interval (0 to 5 seconds, step 0.1)

% Create a figure for the animation
figure;
hold on;
grid on;

% Set axis limits
xlim([min(x), max(x)]);
ylim([-A-0.5, A+0.5]);
xlabel('Position x');
ylabel('Wave Profile eta(x,t)');
title('eta(x,t) = A cos(kx - w t)');

% Plot wave profile for each time step
for t = t_interval
    eta = A * cos(k * x - w * t); % Calculate wave profile
    plot(x, eta, 'b', 'LineWidth', 2);
    
    % Add time stamp as text
    time_str = sprintf('Time: %.2f s', t);
    text(1, A+0.3, time_str, 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
    
    pause(0.1); % Pause to create animation effect
    if t < t_interval(end)
        cla; % Clear axis for next frame
    end
end

hold off;