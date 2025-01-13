clc
d = 10; %depth
T = 3; %Time period
g = 9.81; % gravity
w = (2*pi)/T; %frequency
A = 0.1; %amplitude
%del_k = 0.00004473;
del_w = 0.020944;

% Define the function to solve for k
dispersion_relation = @(k) w^2 - g * k * tanh(k * d);

% Initial guess for k (e.g., shallow water approximation)
k_initial_guess = w^2 / g;

% Solve for k using fzero
k = fzero(dispersion_relation, k_initial_guess);

% Display the result
fprintf('The wave number k is %.4f m^-1\n', k);

% x = linspace(0, 20, 500); % x-axis range
% t_interval = 0:0.1:10;     % Time interval (0 to 10 seconds, step 0.1)
% 
% % Create a figure for the animation
% figure;
% hold on;
% grid on;
% 
% % Set axis limits
% xlim([min(x), max(x)]);
% ylim([-A-0.5, A+0.5]);
% xlabel('Position x');
% ylabel('Wave Profile eta(x,t)');
% title('eta(x,t) = A cos(kx - w t)');
% 
% % Plot wave profile for each time step
% for t = t_interval
%     eta = A * cos(k * x - w * t); % Calculate wave profile
%     plot(x, eta, 'b', 'LineWidth', 2);
% 
%     % Add time stamp as text
%     time_str = sprintf('Time: %.2f s', t);
%     text(1, A+0.3, time_str, 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
% 
%     pause(0.1); % Pause to create animation effect
%     if t < t_interval(end)
%         cla; % Clear axis for next frame
%     end
% end
% 
% hold off;


w1 = w + del_w;
w2 = w - del_w;

disp(w1);
disp(w2);

% Define the function to solve for k
dispersion_relation = @(k) (w1)^2 - g * k * tanh(k * d);

% Initial guess for k (e.g., shallow water approximation)
k1_initial_guess = (w1)^2 / g;

% Solve for k using fzero
k1 = fzero(dispersion_relation, k_initial_guess);

% Display the result
fprintf('The wave number k1 is %.4f m^-1\n', k1);



% Define the function to solve for k
dispersion_relation = @(k) (w2)^2 - g * k * tanh(k * d);

% Initial guess for k (e.g., shallow water approximation)
k2_initial_guess = (w2)^2 / g;

% Solve for k using fzero
k2 = fzero(dispersion_relation, k_initial_guess);

% Display the result
fprintf('The wave number k2 is %.4f m^-1\n', k2);



% % Define spatial and temporal range
% x = linspace(0, 500, 500); % Spatial range
% t_interval = 0:0.1:5;     % Time interval (0 to 5 seconds, step 0.1)
% 
% % Create a figure for the animation
% figure;
% grid on;
% 
% % Set axis limits
% xlim([min(x), max(x)]);
% ylim([-A-A-0.5, A+A+0.5]);
% xlabel('Position x');
% ylabel('Wave Profile \eta(x,t)');
% title('Wave Profiles: \eta_1, \eta_2, and \eta_{new}');
% 
% % Plot wave profiles over time
% for t = t_interval
%     % Calculate wave profiles
%     eta1 = A * cos(k1 * x - w1 * t); % Wave profile of eta1
%     eta2 = A * cos(k2 * x - w2 * t); % Wave profile of eta2
%     eta_new = eta1 + eta2;                % Combined wave profile
% 
%     % Plot the wave profiles
%     plot(x, eta1, 'b', 'LineWidth', 1.5); hold on; % Plot eta1
%     plot(x, eta2, 'r', 'LineWidth', 1.5);          % Plot eta2
%     plot(x, eta_new, 'k', 'LineWidth', 2);         % Plot eta_new
% 
%     % Add legend
%     legend('\eta_1', '\eta_2', '\eta_{new}', 'Location', 'best');
% 
%     % Add time stamp as text
%     time_str = sprintf('Time: %.2f s', t);
%     text(1, A+A+0.3, time_str, 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
% 
%     % Pause for animation effect
%     pause(0.1);
% 
%     % Clear figure for next frame
%     hold off;
% end


% x = linspace(0, 20, 500); % x-axis range
% t_interval = 0:0.1:5;     % Time interval (0 to 10 seconds, step 0.1)
% 
% % Create a figure for the animation
% figure;
% hold on;
% grid on;
% 
% % Set axis limits
% xlim([min(x), max(x)]);
% ylim([-A-0.5, A+0.5]);
% xlabel('Position x');
% ylabel('Wave Profile eta(x,t)');
% title('eta1(x,t) = A cos(k1x - w1t)');
% 
% % Plot wave profile for each time step
% for t = t_interval
%     eta1 = A * cos(k1 * x - w1 * t); % Calculate wave profile
%     subplot(3, 1, 1);
%     plot(x, eta1, 'b', 'LineWidth', 2);
% 
%     % Add time stamp as text
%     time_str = sprintf('Time: %.2f s', t);
%     text(1, A+0.3, time_str, 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
% 
%     pause(0.1); % Pause to create animation effect
%     if t < t_interval(end)
%         cla; % Clear axis for next frame
%     end
% end
% 
% hold off;
% 
% 
% x = linspace(0, 20, 500); % x-axis range
% t_interval = 0:0.1:5;     % Time interval (0 to 10 seconds, step 0.1)
% 
% % Create a figure for the animation
% figure;
% hold on;
% grid on;
% 
% % Set axis limits
% xlim([min(x), max(x)]);
% ylim([-A-0.5, A+0.5]);
% xlabel('Position x');
% ylabel('Wave Profile eta(x,t)');
% title('eta2(x,t) = A cos(k2x - w2t)');
% 
% % Plot wave profile for each time step
% for t = t_interval
%     eta2 = A * cos(k2 * x - w2 * t); % Calculate wave profile
%     subplot(3, 1, 2);
%     plot(x, eta2, 'b', 'LineWidth', 2);
% 
%     % Add time stamp as text
%     time_str = sprintf('Time: %.2f s', t);
%     text(1, A+0.3, time_str, 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');
% 
%     pause(0.1); % Pause to create animation effect
%     if t < t_interval(end)
%         cla; % Clear axis for next frame
%     end
% end
% 
% hold off;

%eta_new = eta1 + eta2;



x = linspace(0, 500, 500); % x-axis range
t_interval = 0:0.1:5;     % Time interval (0 to 10 seconds, step 0.1)

% Create a figure for the animation
figure;
hold on;
grid on;

% Set axis limits
xlim([min(x), max(x)]);
ylim([-A-0.5, A+0.5]);
xlabel('Position x');
ylabel('Wave Profile eta(x,t)');
title('etanew(x,t) = A cos(kx - wt)');

% Plot wave profile for each time step
for t = t_interval
    % Find the crests (maxima) of eta1 and eta2
    [~, idx1] = max(eta1); % Index of crest for eta1
    [~, idx2] = max(eta2); % Index of crest for eta2
    subplot(3,1,1);
    plot(x, eta1, 'b', 'LineWidth', 2);
    eta1 = A * cos(k1 * x - w1 * t);
    subplot(3,1,2);
    plot(x, eta2, 'b', 'LineWidth', 2);
    eta2 = A * cos(k2 * x - w2 * t); % Calculate wave profile
    eta_new = eta1 + eta2;
    subplot(3, 1, 3);
    plot(x, eta_new, 'b', 'LineWidth', 2);

    % Add time stamp as text
    time_str = sprintf('Time: %.2f s', t);
    text(1, A+0.3, time_str, 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'r');

    pause(0.1); % Pause to create animation effect
    if t < t_interval(end)
        cla; % Clear axis for next frame
    end
end

hold off;


