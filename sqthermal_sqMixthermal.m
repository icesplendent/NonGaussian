clc; clear all;close all;
%%  squeezed thermal state
% parameters
rdB = 5; %sq dB
r = rdB/20*log(10); %squeezing factor
phi = 0; %sq angle  *phase=0¬°amplitude-sq, phase=pi¬°phase-sq
alpha = 0; %coherent
phase = linspace(0,2*pi,200);
N = 100; % should be larger than alpha&squeezing factor
Bol = 0.5; %(exp(hbar*w/k_b*T)-1)^(-1)
% Define operators
a0 = diag(sqrt(1:N),1); %annihilation  operator
ad0 = a0'; %creation operator
Dis = expm(alpha*a0'-conj(alpha)*a0); %displacement operator
S = expm(0.5*(r*exp(-1i*phi)*a0^2-r*exp(1i*phi)*ad0^2)); %squeezing operator  2*

% Prepare state
thermal = zeros(N+1,N+1);
for nn = 0:N
    thermal(nn+1,nn+1) = Bol^nn/((1+Bol)^(nn+1));
end
state = Dis*S*thermal*S'*Dis';
figure;bar3(state/trace(state))
%% squeezed state + thermal state
vacuum = zeros(N+1,N+1);
vacuum(1,1) = 1;

thermal = zeros(N+1,N+1);
for nn = 0:N
    thermal(nn+1,nn+1) = Bol^nn/((1+Bol)^(nn+1));
end
state_sq = Dis*S*vacuum*S'*Dis';
state_mix = state_sq + thermal/trace(thermal);
figure;bar3(state_mix/trace(state_mix))