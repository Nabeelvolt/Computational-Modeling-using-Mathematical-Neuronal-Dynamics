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
T = 5;                 % Total simulation time (s)
dt = 0.1e-3;           % Time step (s)
t = 0:dt:T;            % Time vector

% Applied current range
I_min = 0;             % Minimum current (nA)
I_max = 0.5;             % Maximum current (nA)
n_currents = 20;       % Number of current levels
I_range = linspace(I_min, I_max, n_currents) * 1e-9; % Current range in amperes

% Variables to store results
f_initial = zeros(1, n_currents);  % Inverse of initial ISI (Hz)
f_steady = zeros(1, n_currents);   % Inverse of steady-state ISI (Hz)

% Loop over applied current values
for j = 1:n_currents
    % Initialize variables
    V = E_L * ones(1, length(t));    % Membrane potential (V)
    I_SRA = zeros(1, length(t));     % Adaptation current (A)
    spike_times = [];                % Spike times

    % Simulate for current level I_range(j)
    for i = 2:length(t)
        % Update adaptation current
        dI_SRA = (a * (V(i-1) - E_L) - I_SRA(i-1)) / tau_SRA;
        I_SRA(i) = I_SRA(i-1) + dI_SRA * dt;

        % Update membrane potential using Euler's method
        dV = (G_L * (E_L - V(i-1) + Delta_th * exp((V(i-1) - V_th) / Delta_th)) - I_SRA(i) + I_range(j)) / C_m;
        V(i) = V(i-1) + dV * dt;

        % Check for spike
        if V(i) >= V_th
            V(i) = V_reset;          % Reset membrane potential
            I_SRA(i) = I_SRA(i) + b; % Increment adaptation current
            spike_times = [spike_times, t(i)]; % Record spike time
        end
    end

    % Calculate ISIs
    ISIs = diff(spike_times); % Time differences between spikes

    % Initial ISI (first spike interval)
    if length(ISIs) >= 1
        f_initial(j) = 1 / ISIs(1); % Inverse of the first ISI (Hz)
    else
        f_initial(j) = 0; % No spikes occurred
    end

    % Steady-state ISI (average ISI in the last second)
    steady_spikes = spike_times(spike_times > T - 1); % Spikes in the last second
    if length(steady_spikes) >= 2
        steady_ISIs = diff(steady_spikes); % Steady-state ISIs
        f_steady(j) = 1 / mean(steady_ISIs); % Inverse of the steady-state ISI (Hz)
    else
        f_steady(j) = 0; % No steady-state spikes
    end
end

% Plotting the f-I curve
figure;
plot(I_range * 1e9, f_steady, '-o', 'LineWidth', 2); % Steady-state
hold on;
plot(I_range * 1e9, f_initial, '-x', 'LineWidth', 2); % Initial
xlabel('Applied Current (nA)', 'FontSize', 14);
ylabel('Firing Rate (Hz)', 'FontSize', 14);
title('f-I Curve (AELIF Model)', 'FontSize', 16);
legend('Steady-State', 'Initial', 'FontSize', 12);
grid on;
