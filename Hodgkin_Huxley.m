% Hodgkin-Huxley Neuron Model Simulation
clear all;

% Parameters
C_m = 100e-12;    % Membrane capacitance (F/cm^2)
g_Na = 12e-6;     % Maximum sodium conductance (S/cm^2)
g_K = 3.6e-6;     % Maximum potassium conductance (S/cm^2)
g_L = 30e-9;      % Leak conductance (S/cm^2)
E_Na = 45e-3;     % Sodium reversal potential (V)
E_K = -82e-3;     % Potassium reversal potential (V)
E_L = -60e-3;     % Leak reversal potential (V)
V_rest = -70.2e-3;% Resting membrane potential (V)

% Time parameters
dt = 0.01e-3;  % Time step (s)
T = 0.35;      % Total simulation time (s)
time = 0:dt:T; % Time vector

% External current (stimulation)
I_ext = zeros(size(time));
I_ext(ceil(0.05/dt):ceil(0.06/dt)) = 0;  % No external current

% Initialize variables
V = V_rest * ones(size(time));  % Membrane potential (V)
m = 0.0;   % Sodium activation gating variable
h = 0.0;   % Sodium inactivation gating variable
n = 0.0;   % Potassium activation gating variable

% Storage for gating variables over time
m_values = zeros(size(time));
h_values = zeros(size(time));
n_values = zeros(size(time));
V_values = zeros(size(time));

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
    m = m + dt * (alpha_m(V(t-1)) * (1 - m) - beta_m(V(t-1)) * m);
    h = h + dt * (alpha_h(V(t-1)) * (1 - h) - beta_h(V(t-1)) * h);
    n = n + dt * (alpha_n(V(t-1)) * (1 - n) - beta_n(V(t-1)) * n);
    
    % Compute conductances for sodium and potassium
    g_Na_t = g_Na * (m^3) * h;
    g_K_t = g_K * (n^4);
    
    % Compute ionic currents
    I_Na = g_Na_t * (V(t-1) - E_Na);
    I_K = g_K_t * (V(t-1) - E_K);
    I_L = g_L * (V(t-1) - E_L);
    
    % Update membrane potential using Euler's method
    V(t) = V(t-1) + dt * (I_ext(t) - (I_Na + I_K + I_L)) / C_m;
    
    % Store values for plotting
    m_values(t) = m;
    h_values(t) = h;
    n_values(t) = n;
    V_values(t) = V(t);
end

% Plot membrane potential over time
figure;
subplot(2,1,1);
plot(time * 1e3, V_values * 1e3, 'LineWidth', 2);
title('Hodgkin-Huxley Model: Membrane Potential');
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
grid on;

% Plot gating variables over time
subplot(2,1,2);
plot(time * 1e3, m_values, 'r', 'LineWidth', 1.5); hold on;
plot(time * 1e3, h_values, 'g', 'LineWidth', 1.5);
plot(time * 1e3, n_values, 'b', 'LineWidth', 1.5);
legend('m (Na+ activation)', 'h (Na+ inactivation)', 'n (K+ activation)');
title('Gating Variables Over Time');
xlabel('Time (ms)');
ylabel('Gating Variable');
grid on;
