% Parameters
E_L = -75e-3;           % Resting potential (mV)
V_th = -50e-3;          % Threshold potential (mV)
V_reset = -80e-3;       % Reset potential (mV)
R_m = 100e6;           % Membrane resistance (MOhm)
C_m = 100e-12;           % Membrane capacitance (pF)
E_k = -80e-3;           % Adaptation reversal potential (mV)
tau_SRA = 200e-3;       % Adaptation time constant (ms)
Delta_G_SRA = 1e-9;  % Spike-triggered adaptation increment (uS)
I_app = 5e-10;         % Applied current (nA)
dt = 0.1e-3;            % Time step (ms)
T = 1.5;            % Total simulation time (ms)

% Time vector
time = 0:dt:T;
n_steps = length(time);

% Variables
V = E_L * ones(1, n_steps);      % Membrane potential (mV)
G_SRA = zeros(1, n_steps);       % Adaptation conductance (uS)
spike_train = zeros(1, n_steps); % Spike train

for t = 2:n_steps
    % Update G_SRA (adaptation conductance)
    dG_SRA = -G_SRA(t-1) / tau_SRA;
    G_SRA(t) = G_SRA(t-1) + dG_SRA * dt;

    % Update membrane potential using Euler's method
    dV = ((E_L - V(t-1)) / R_m + G_SRA(t) * (E_k - V(t-1)) + I_app) * (dt / C_m);
    V(t) = V(t-1) + dV;

    % Check for spike
    if V(t) > V_th
        V(t) = V_reset; % Reset potential
        G_SRA(t) = G_SRA(t) + Delta_G_SRA; % Increment adaptation conductance
        spike_train(t) = 1; % Record spike
    end
end

% Subplot 1: Membrane potential
subplot(3, 1, 1);
plot(time, V, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Membrane Potential (V)');
title('LIF Model with Adaptation: Membrane Potential');
grid on;

% Subplot 2: Adaptation conductance
subplot(3, 1, 2);
plot(time, G_SRA, 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Adaptation Conductance (S)');
title('LIF Model with Adaptation: Adaptation Conductance');
grid on;

% Subplot 3: Spike train
subplot(3, 1, 3);
stem(time, spike_train, 'Marker', 'none');
xlabel('Time (s)');
ylabel('Spike');
title('LIF Model with Adaptation: Spike Train');
grid on;