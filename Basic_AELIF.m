% Parameters for the AELIF model (in SI units)
E_L = -75e-3;          % Resting membrane potential (V)
V_th = -50e-3;         % Threshold potential (V)
V_reset = -80e-3;      % Reset potential after spike (V)
Delta_th = 2e-3;       % Exponential slope factor (V)
G_L = 10e-9;           % Leak conductance (S)
C_m = 100e-12;         % Membrane capacitance (F)
a = 2e-9;              % Subthreshold adaptation conductance (S)
b = 0.02e-9;           % Spike-triggered adaptation increment (A)
tau_SRA = 200e-3;      % Adaptation time constant (s)
T = 0.5;               % Total simulation time (s)
dt = 0.1e-3;           % Time step (s)
t = 0:dt:T;            % Time vector

% Initial conditions
V = E_L * ones(1, length(t));    % Membrane potential (V)
I_SRA = zeros(1, length(t));     % Adaptation current (A)
spike_train = zeros(1, length(t)); % Spike train

% Loop to simulate the AELIF model
for i = 2:length(t)
    % Update adaptation current
    dI_SRA = (a * (V(i-1) - E_L) - I_SRA(i-1)) / tau_SRA;
    I_SRA(i) = I_SRA(i-1) + dI_SRA * dt;

    % Update membrane potential using Euler's method
    dV = (G_L * (E_L - V(i-1) + Delta_th * exp((V(i-1) - V_th) / Delta_th)) - I_SRA(i)) / C_m;
    V(i) = V(i-1) + dV * dt;

    % Check for spike
    if V(i) >= V_th
        V(i) = V_reset;          % Reset membrane potential
        I_SRA(i) = I_SRA(i) + b; % Increment adaptation current
        spike_train(i) = 1;      % Record spike
    end
end

% Plotting the results
figure;

% Membrane Potential
subplot(2, 1, 1);
plot(t, V * 1e3, 'LineWidth', 2); % Convert back to mV for display
xlabel('Time (s)', 'FontSize', 14);
ylabel('Membrane Potential (mV)', 'FontSize', 14);
title('Membrane Potential (AELIF Model)', 'FontSize', 16);
grid on;

% Adaptation Current
subplot(2, 1, 2);
plot(t, I_SRA * 1e9, 'LineWidth', 2); % Convert back to nA for display
xlabel('Time (s)', 'FontSize', 14);
ylabel('Adaptation Current (nA)', 'FontSize', 14);
title('Adaptation Current', 'FontSize', 16);
grid on;
