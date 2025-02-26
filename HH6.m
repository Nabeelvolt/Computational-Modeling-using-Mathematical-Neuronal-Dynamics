% Hodgkin-Huxley Neuron Model Simulation for Question 6
clear all;

% Parameters
C_m = 100e-12;    % Membrane capacitance (F/cm^2)
g_Na = 12e-6;     % Maximum sodium conductance (S/cm^2)
g_K = 3.6e-6;     % Maximum potassium conductance (S/cm^2)
g_L = 30e-9;      % Leak conductance (S/cm^2)
E_Na = 45e-3;     % Sodium reversal potential (V)
E_K = -82e-3;     % Potassium reversal potential (V)
E_L = -60e-3;     % Leak reversal potential (V)
V_rest = -65e-3;  % Resting membrane potential (V)

% Initial conditions
V = V_rest;       % Membrane potential (V)
m = 0;            % Sodium activation gating variable
h = 0;            % Sodium inactivation gating variable
n = 0;            % Potassium activation gating variable

% Time parameters
dt = 0.01e-3;  % Time step (s)
T = 0.5;       % Total simulation time (s)
time = 0:dt:T; % Time vector

% Baseline applied current (0.7 nA)
baseline_current = 0.7e-9; % Baseline current (A)
I_ext = baseline_current * ones(size(time));

% Add a single excitatory pulse of 1 nA at 100 ms
pulse_amplitude = 1e-9;  % Excitatory pulse amplitude (A)
pulse_duration = 5e-3;   % Pulse duration (s)
pulse_time = 0.1;        % Pulse start time (s)

% Define pulse start and end indices
pulse_start_idx = max(1, ceil(pulse_time / dt));
pulse_end_idx = min(length(time), pulse_start_idx + ceil(pulse_duration / dt) - 1);

% Apply excitatory pulse
I_ext(pulse_start_idx:pulse_end_idx) = I_ext(pulse_start_idx:pulse_end_idx) + pulse_amplitude;

% Initialize variables for simulation
V_values = zeros(size(time));
m_values = zeros(size(time));
h_values = zeros(size(time));
n_values = zeros(size(time));

% Set initial conditions
V_values(1) = V;
m_values(1) = m;
h_values(1) = h;
n_values(1) = n;

% Define rate functions for gating variables
alpha_m = @(V) (1e5 * (V + 0.045)) ./ (exp(100 * (V + 0.045)) - 1);
beta_m = @(V) 4e3 * exp((V + 0.070) / 0.018);
alpha_h = @(V) 70 * exp(50 * (V + 0.070));
beta_h = @(V) 1e3 ./ (1 + exp(100 * (V + 0.040)));
alpha_n = @(V) (1e4 * (V + 0.060)) ./ (exp(100 * (V + 0.060)) - 1);
beta_n = @(V) 125 * exp((V + 0.070) / 0.08);

% Simulation loop
for t = 2:length(time)
    % Update gating variables using Euler method
    m = m + dt * (alpha_m(V_values(t-1)) * (1 - m) - beta_m(V_values(t-1)) * m);
    h = h + dt * (alpha_h(V_values(t-1)) * (1 - h) - beta_h(V_values(t-1)) * h);
    n = n + dt * (alpha_n(V_values(t-1)) * (1 - n) - beta_n(V_values(t-1)) * n);
    
    % Compute conductances for sodium and potassium
    g_Na_t = g_Na * (m^3) * h;
    g_K_t = g_K * (n^4);
    
    % Compute ionic currents
    I_Na = g_Na_t * (V_values(t-1) - E_Na);
    I_K = g_K_t * (V_values(t-1) - E_K);
    I_L = g_L * (V_values(t-1) - E_L);
    
    % Update membrane potential using Euler's method
    V_values(t) = V_values(t-1) + dt * (I_ext(t) - (I_Na + I_K + I_L)) / C_m;
    
    % Store gating variable values
    m_values(t) = m;
    h_values(t) = h;
    n_values(t) = n;
end

% Plot applied current and membrane potential
figure;

% Plot applied current
subplot(2, 1, 1);
plot(time * 1e3, I_ext * 1e9, 'LineWidth', 2); % Convert to nA for plotting
title('Applied Current with Excitatory Pulse (Baseline: 0.7 nA)');
xlabel('Time (ms)');
ylabel('Current (nA)');
grid on;

% Plot membrane potential
subplot(2, 1, 2);
plot(time * 1e3, V_values * 1e3, 'LineWidth', 2); % Convert to mV for plotting
title('Membrane Potential (V_m)');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
grid on;
