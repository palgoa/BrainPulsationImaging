function [S1,S2,X,V]  = EPG_nbSSFP_brainpulsation_sinus(sim)
%
% function [S1,S2,X,V]  = EPG_ssfp_brainpuls(sim)
%
% Simple EPG simulation on the effect of sinusoidal motion on the signal in
% non-balanced SSFP sequence. TE is set to zero.
% 
% All units SI
% inputs: (no error-handling implemented).
% sim must contain following item:
% sim.RR - period of motion [sec].
% sim.Xx - amplitude of motion [m].
% sim.T1 - Longitudinal Relaxation Time [sec]
% sim.T2 - transverse relaxation time [sec]
% sim.M0 - Gradient spoiler moment at the end of the TR-interval
% (immediately before the next RF-pulse).
% sim.M1 - First order Gradient spoiler moment at the end of the TR-interval
% (immediately before the next RF-pulse).
% sim.FA - RF-pulse flip angle [deg] 
% sim.np -duration of simulation (number of TR-periods).
%
% Pål Erik Goa NTNU 2025
% pal.e.goa@ntnu.no
%
gamma = 267.5e6; % rad/T/s

%Period of displacement and velocity oscillations.
Tv = sim.RR; 

% Displacement amplitudes
Xx = sim.Xx; 

%Velocity amplitudes (from time derivative of Xx*sin(2*pi*t/Tv);
Vx = Xx*(2*pi/Tv); 

% Material parameters:
T1 = sim.T1;
T2 = sim.T2;

% Sequence Parameters:
TR = sim.TR;
alpha = sim.FA*pi/180; %rad
M0  = sim.M0;
M1  = sim.M1;

% Duration of simulation in terms of number of TRs:
np = sim.np; 
% Limits the number of states kept in the simulation:
maxstates = 200;

% Terms used later:
E1 = exp(-TR/T1);
E2 = exp(-TR/T2);
a = (cos(alpha/2))^2;
b = (sin(alpha/2))^2;
c = sin(alpha);
d = cos(alpha);

phi = 0; % No RF-spoiling
e = exp(1i*2*phi);
f = exp(-1i*2*phi);
g = exp(1i*phi);
h = exp(-1i*phi);

%Initialzation of variables:
F = zeros(2*np-1,1);
pF = F;
Z = zeros(np-1,1);
pZ = Z;
S1 = Z;
S2 = Z;
X=Z;
V = Z;

Z(1)  = 1; % initial state

% for F and pF states: Zero order state is in index np;
% For z-states: Zero order state is in index 1;

for k = 1:np-1 % for all pulses
    X(k) = Xx*sin(2*pi*k*TR/Tv);
    V(k) = Vx*cos(2*pi*k*TR/Tv); 
   
    % Effect of RF-pulse:
    for n = 0:min(k-1,maxstates) % for all states created so far
        npos =np+n; % index for positive n;
        nneg =np-n; % index for negative n;
        nz = n+1;
        pF(npos) = a*F(npos) + e*b*conj(F(nneg)) - 1i*g*c*Z(nz);
        pF(nneg) = conj(f*b*F(npos) + a*conj(F(nneg)) + 1i*h*c*Z(nz));
        pZ(nz) = -1i/2*h*c*F(npos) + 1i/2*g*c*conj(F(nneg)) + d*Z(nz);
    end
    % Calculating global phase-term at end of TR due to motion:
    TRphase(k) = gamma*(M0*X(k)+ M1*V(k));

    % Propagation through TR:
    for n = 0:min(k-1,maxstates) % for all states
        npos =np+n; % index for positive n;
        nneg =np-n; % index for negative n;
        F(npos+1) = pF(npos) * E2 *exp(1i*TRphase(k)); % Non-balanced
        F(nneg+1) = pF(nneg) * E2 *exp(1i*TRphase(k)); % Non-balanced
        if n > 0
            Z(n+1) = E1*pZ(n+1);
        else
            Z(1) = E1*pZ(1) + 1*(1-E1);
        end
    end
    S1(k) = pF(np);
    S2(k) = F(np);
end




