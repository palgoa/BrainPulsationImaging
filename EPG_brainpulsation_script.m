% Script to run EPG-simulations for non-balanced SSFP with sinusoidal
% displacement.
%
% Supplementary code for doi:10.1002/mrm.70156
%
% Pål Erik Goa 2025
%
% 

%% Material parameters at 7T
% Uncomment the tissue type you want to simulate:
%
% CSF:
sim.T1 = 4.4; %Rooney 2007 /MArques 2018
sim.T2 = 1.0; % Spijkerman 2018
%
% WM Tissue
% sim.T1 = 1.2; %seconds Wright 2008
% sim.T2 = 0.055; % Bartha 2002/ Marques 2018
%
%GM Tissue
% sim.T1 = 2.0; %seconds Wright 2008
% sim.T2 = 0.05; % Marques 2018

%% Motion Parameters:
%
% Period of sinusoidal displacement:
sim.RR = 1.0; % [sec]
% Amplitude of sinusoidal displacement:
sim.Xx = [1e-3]; % [m]

%% Sequence parameters
sim.TR = 10e-3; % [sec]
sim.FA = 20;    % [deg] 
% Gradient spoiler moments (values are correct for the 3D-EPI sequence used in
% brainpulsation imaging paper for 2pi zero order moment
sim.M0 = 7.339e-6; % [Ts/m]
sim.M1 = 68.7e-9;  % [Ts^2/m]

%% EPG Parameters:
sim.np = round(5*sim.T1/sim.TR); % At least 5 times T1 to reach "steady state"

%% Run simulation
[S1 S2 X V] = EPG_nbSSFP_brainpulsation_sinus(sim);

%% Plot results 

% Last period of the oscillation:
OnePeriod = round(sim.RR./sim.TR);

figure, subplot(3,1,1), plot([0:OnePeriod].*sim.TR,1e3.*X(end-OnePeriod:end),LineWidth=2);
title('Displacement')
xlabel('Time [sec]')
ylabel('Displacement [mm]')

subplot(3,1,2), plot([0:OnePeriod].*sim.TR,abs(S1(end-OnePeriod:end)),LineWidth=2);
title('S1 Magnitude')
xlabel('Time [sec]')
ylabel('nbSSFP Magnitude[a.u.]')

S1RelPhase = unwrap(angle(S1))- mean(unwrap(angle(S1)));
subplot(3,1,3), plot([0:OnePeriod].*sim.TR,S1RelPhase(end-OnePeriod:end),LineWidth=2)
title('S1 Relative Phase')
xlabel('Time [sec]')
ylabel('nbSSFP Phase [rad]')
