% Parameters of the LIF model with adaptation (in SI units)
E_L = -75e-3;          % Resting potential (V)
V_th = -50e-3;         % Threshold potential (V)
V_reset = -80e-3;      % Reset potential after spike (V)
R_m = 100e6;           % Membrane resistance (Ohms)
C_m = 100e-12;         % Membrane capacitance (Farads)
E_k = -80e-3;          % Adaptation reversal potential (V)
tau_SRA = 200e-3;      % Adaptation time constant (s)
Delta_G_SRA = 1e-9;    % Increment in adaptation conductance (Siemens)
T = 5;                 % Total simulation time (s)
dt = 0.1e-3;           % Time step (s)
t = 0:dt:T;            % Time vector

% Applied current range
I_min = 0;             % Minimum current (nA)
I_max = 0.5;           % Maximum current (nA)
n_currents = 20;       % Number of current levels
I_range = linspace(I_min, I_max, n_currents) * 1e-9; % Current range in amperes

% Variables to store results
f_initial = zeros(1, n_currents);  % Inverse of initial ISI
f_steady = zeros(1, n_currents);   % Inverse of steady-state ISI

% Loop over applied current values
for j = 1:n_currents
    % Initialize variables
    V = E_L * ones(1, length(t));    % Membrane potential (V)
    G_SRA = zeros(1, length(t));     % Adaptation conductance (S)
    spike_times = [];                % Spike times

    % Simulate for current level I_range(j)
    for i = 2:length(t)
        % Update adaptation conductance
        dG_SRA = -G_SRA(i-1) / tau_SRA;
        G_SRA(i) = G_SRA(i-1) + dG_SRA * dt;

        % Update membrane potential using Euler's method
        dV = ((E_L - V(i-1)) / R_m + G_SRA(i) * (E_k - V(i-1)) + I_range(j)) * (dt / C_m);
        V(i) = V(i-1) + dV;

        % Check for spike
        if V(i) >= V_th
            V(i) = V_reset;  % Reset the potential
            G_SRA(i) = G_SRA(i) + Delta_G_SRA; % Increase adaptation conductance
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
title('f-I Curve (LIF Model with Adaptation)', 'FontSize', 16);
legend('Steady-State', 'Initial', 'FontSize', 12);
grid on;