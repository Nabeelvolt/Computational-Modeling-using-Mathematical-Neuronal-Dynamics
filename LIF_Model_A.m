% Parameters of the LIF model with adaptation (in SI units)
E_L = -75e-3;          % Resting potential (V)
V_th = -50e-3;         % Threshold potential (V)
V_reset = -80e-3;      % Reset potential after spike (V)
R_m = 100e6;           % Membrane resistance (Ohms)
C_m = 100e-12;         % Membrane capacitance (Farads)
E_k = -80e-3;          % Adaptation reversal potential (V)
tau_SRA = 200e-3;      % Adaptation time constant (s)
Delta_G_SRA = 1e-9;    % Increment in adaptation conductance (Siemens)
T = 1.5;               % Total simulation time (s)
dt = 0.1e-3;           % Time step (s)
t = 0:dt:T;            % Time vector

% Input current (in amperes)
I_app = zeros(1, length(t));      % Initialize input current
I_app(0.5/dt:1.0/dt) = 500e-12;  % Apply 500 pA current from 0.5s to 1.0s

% Variables
V = E_L * ones(1, length(t));    % Membrane potential (V)
G_SRA = zeros(1, length(t));     % Adaptation conductance (S)
spike_train = zeros(1, length(t)); % Spike train

% Simulating the LIF neuron with adaptation
for i = 2:length(t)
    % Update adaptation conductance
    dG_SRA = -G_SRA(i-1) / tau_SRA;
    G_SRA(i) = G_SRA(i-1) + dG_SRA * dt;

    % Update membrane potential using Euler's method
    dV = ((E_L - V(i-1)) / R_m + G_SRA(i) * (E_k - V(i-1)) + I_app(i)) * (dt / C_m);
    V(i) = V(i-1) + dV;

    % Check for spike
    if V(i) >= V_th
        V(i) = V_reset;            % Reset the potential
        spike_train(i) = 1;        % Record a spike
        G_SRA(i) = G_SRA(i) + Delta_G_SRA; % Increase adaptation conductance
    end
end

% Plotting the results
figure;

% Applied Current
subplot(3, 1, 1);
plot(t, I_app * 1e12, 'LineWidth', 2); % Convert back to pA for display
xlabel('Time (ms)', 'FontSize', 14);
ylabel('Current (pA)', 'FontSize', 14);
title('Applied Current', 'FontSize', 16);
grid on;

% Membrane Potential
subplot(3, 1, 2);
plot(t, V * 1e3, 'LineWidth', 2); % Convert back to mV for display
xlabel('Time (ms)', 'FontSize', 14);
ylabel('Membrane Potential (mV)', 'FontSize', 14);
title('Membrane Potential with Adaptation', 'FontSize', 16);
ylim([-85 -45]);
grid on;

% Adaptation Conductance
subplot(3, 1, 3);
plot(t, G_SRA * 1e9, 'LineWidth', 2); % Convert back to nS for display
xlabel('Time (ms)', 'FontSize', 14);
ylabel('Adaptation Conductance (nS)', 'FontSize', 14);
title('Adaptation Conductance', 'FontSize', 16);
grid on;
