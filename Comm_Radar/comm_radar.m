%% For comparing with UE - AP data
clc;
clear;
%function main

%% Communication parameters
% AP parameters
APAngR = [-40];
ref= [ 20 ]; %
ntarg = length(APAngR); %number of targets

M = 4;                              % n APs
NAp = 200;                           % nr antennas Acess Point (communication)
Q = NAp;                            %nr antennas for receiving AP (radar)
P = NAp;                              %nr antennas for transmitting AP (radar)

% UE  parameters
U = 4;                              % nr UEs
NUe = 1;                           % nr antennas User Equipment (communication)
NrxRF = NUe;                        % nr RF chains is equal to the number of antennas

% channel parameters

Ncl = 5;                            % nr clusters
Nray = 3;                           % nr rays
AngleSpread = 8;                    % angle spread

% OFDM parameters
K = 64;                            % nr of subcarriers
MODULATION_ORDER = 4;               % QPSK

% simulation parameters
Nchannels = 10;                    % number of channel realizations
Nsym = 4;                          % number of noise realizations
EbN0 = -28:4:4;                   % energy per bit to noise power spectral density ratio 
%EbN0=100;
nvar = 10.^(-EbN0./10)/(log2(MODULATION_ORDER));   % noise variance

%% Radar Parameters
c0 = 3e8;
fc = 28e9;  % carrier frequency (originalmente 24 mas e para ser consistente com o que e utilizado aqui)
lambda = 3e8 / fc;

B = 93.1e6; % signal bandwidth
Tsymbol = K / B; % Useful symbol duration
Tcp = (1 / 4) * Tsymbol;
To = Tsymbol + Tcp;

spacing = 1 / Tsymbol;
d = lambda / 2;         % Element spacing
kw = 2 * pi / lambda;    % Wave number
elemR = (0:NAp-1)';

%% allocate memory and construct objects
ber = nan(Nchannels, Nsym, length(EbN0));
channel = CWideband_mmWave_Channel(NUe, NAp, K, K/4, Ncl, Nray, 10, 10, AngleSpread); %% construct channel object (DL)
DFT_S_OFDM = CDFT_S_OFDM(K, U, 1, MODULATION_ORDER); %% construct DFT-S-OFDM object //sc-fdma (DL)

    for h = 1:ntarg
        for m = 1:M
            APAngT(m,h) = genAngle();
        end   
    end

Fm = zeros(NAp, NAp); % 3D matrix to store each Fm result
AtTotal = zeros(NAp, M); % 2D matrix to store each At result

% calcular o canal 
for m = 1: M
    for h = 1 : ntarg
        [Ar, At] = genRadarChannel(c0, fc, NAp, NAp, APAngT(m,h), APAngR(h) );
        %adicionar ganho de refletividade
        Fm(:,:,m,h) = (Ar*At')*10^(ref(h)./20); 
        AtTotal(:,m,h) = At;
    end 
end 

%% open new figure to plot the results
%% start timer 
tic

% Initialize periodogram accumulation
angles = -90:0.1:90; % Angle grid
P_accum = zeros(size(angles));

%% for each channel realization
for ch = 1:Nchannels

    % generate channel for all UT and AP pairs
    H0 = generateChannels(channel, U, M); %(DL)
    %% for each noise realization
    for n = 1:Nsym
        % generate noise with unit variance
        noise = sqrt(1/2)*complex(randn(NUe, K), randn(NUe, K)); %(DL)
        %% sweep over the EbN0 range
        for p = 1:length(EbN0)
            %% Access Points
            b = DFT_S_OFDM.genRandBits();                                       % Generate random bits
            s = DFT_S_OFDM.mod(b);                                              % DFT_S_OFDM modulation

            H = permute(H0,[4 2 3 5 1]);

            Wd = computeDigitalPrecoder(H, nvar(p), M);

            y1 = zeros(K, U);
            y = zeros(K, U); 
            
            for k = 1:K
                for u = 1:U
                    for m = 1:M
                        y1(k, u) = y1(k, u) + H0(:, :, k, u, m) * Wd(:, :, k, m) * s(k, :).' ;   
                    end
                    y(k, u) = y1(k, u) + sqrt(nvar(p)) * noise(:, k);
                end
            end
                     
            % br = DFT_S_OFDM.demod(y);
            % [~, ber(ch, n, p)] = biterr(b, br);
            
            

        end
    end
    sim_time = toc; % read timer

    %% display information about simulation
    if(mod(ch, 10) == 0)
        avg_sim_time = sim_time/ch;
        Nchannel_left = Nchannels - ch;
        fprintf('Simulation time left (estimate): %0.1f min. \n', avg_sim_time*Nchannel_left/60)
    end
end

% Fixed noise variance
fixed_noise_var = 10^(30/10);
fixed_noise = sqrt(fixed_noise_var/2) * (randn(NAp, 1) + 1j * randn(NAp, 1));

for h = 1 : ntarg
    for m = 1:M
        xm = AtTotal(:,m,h); 
        yr(:,m,h) = Fm(:,:,m,h)*xm ;  
    end
end


% Summing received signals
RxSignal = sum(sum(yr, 3), 2) + fixed_noise;

N = length(RxSignal); % Length of RxSignal
angles = 0:180; % Angle grid for estimation
P = zeros(size(angles));

S = fft(RxSignal, N);


% Precompute steering vectors for all angles
steering_vectors = zeros(NAp, length(angles));
for idx = 1:length(angles)
    angle = angles(idx);
    steering_vectors(:, idx) = exp(1j * elemR * kw * d * sind(angle-90)) / sqrt(NAp);
end

% FFT of all steering vectors
W = fft(steering_vectors, N);

% Compute periodogram using cross-power spectral density approach
for idx = 1:length(angles)
    P(idx) = ((1/ N) * abs(W(:, idx)' * S)^2)/NAp;
end

% Normalize periodogram
%P = P / max(P);
xaxis = -90:90;

% Plot periodogram
figure;
plot(xaxis, 10*log10(P));
xlabel('Angle (degrees)');
ylabel('Normalized Power');
title('Periodogram for Angle Estimation (FFT)');
grid on;


%% generate equivalent channel, including precoder, for all UT and AP pairs
function H0 = generateChannels(channel, M, U)
%% Channel parameters
np = 4.1;   % path-Loss exponents
dp_s = 7.6; % standard deviation of shadowing factor
f = 28;     % GHz
rAP = 500;  % AP coverage radius

%% allocate memory for variables
H0 = zeros(channel.Nr, channel.Nt, channel.Nc, M, U);
Arx0 = zeros(channel.Nr, channel.Ncl * channel.Nray, M, U);
pl = zeros(U, M);
muTx = zeros(U, channel.Ncl, M);

%% generate users and access points positions
[xa, ya, xu, yu] = genNodesPositions(rAP, M, U);

%% Average path loss for a scenario with a user at cell border being served by an AP at cell center
averageShadowing = exp((log(10)/10 * dp_s)^2 / 2);  % theoretical value of average shadowing
pathLossWithoutShadowing = pLoss(f, rAP, np, 0);   % path loss not considering shadowing
C = pathLossWithoutShadowing * averageShadowing;  % normalization constant

% for each AP
for m = 1:M
    % for each UT
    for u = 1:U
        % generate channel without path loss
        [H0(:, :, :, m, u), ~, Arx0(:, :, m, u), ~, ~, muTx(u, :, m)] = channel.genRandomInstance;

        % distance (D) between UTs and APs:
        D = mdist([xu(u), yu(u)], [xa(m), ya(m)]);

        % Path Loss
        pl(u, m) = pLoss(f, D, np, dp_s) / C;

        H0(:, :, :, m, u) = sqrt(pl(u, m)) * H0(:, :, :, m, u);
    end
end
end

function Wd = computeDigitalPrecoder(H, nvar, M)
    % Get parameters
    NAp = size(H, 2);  % Number of receivers
    U = size(H, 1);    % Number of users
    K = size(H, 3);    % Number of subcarriers

    % Allocate memory for Gd
    Wd = zeros(NAp, U, K, M);

    % Compute Gd without Ga
    for m = 1:M
        for k = 1:K
            % Compute received signal correlation
            R = H(:, :, k, m) * H(:, :, k, m)' + nvar * eye(U);
            % Compute optimum digital part
            Wd(:, :, k, m) = H(:, :, k, m)' / R;
            powers = diag(Wd(:, :, k, m)' * Wd(:, :, k, m));
            Wd(:, :, k, m) = Wd(:, :, k, m) / sqrt(M * diag(powers));
        end
    end
end

%% Path loss
function PL = pLoss(f, d, E_NLOS, f_NLOS)
waveLen = 0.3 / f;

G_TX = 15;      % 15 dBi
G_Rx = 24.5;    % 24.5 dBi
Gain = G_TX + G_Rx;
Gain_L = 10^(Gain / 10);

d0 = 1;
Beta0 = 10 * log10((4 * pi * d0) / waveLen)^2;

%E_NLOS = 4.1;  % NLOS Path-loss exponents in dB
%f_NLOS = 7.6;  % NLOS Standard deviation of shadowing factor in dB

A_NLOS = f_NLOS * randn(); % A_NLOS= 10.^(f_NLOS*randn()/10);  % linear

Beta = Beta0 + 10 * E_NLOS * log10(d / d0) + A_NLOS;
Beta_L = 10^(Beta / 10);

%Ga_Beta = Gain - Beta;
%Ga_Beta_L1 = 10^(Ga_Beta / 10);

PL = Gain_L / Beta_L;
end

% Generate nodes positions
function [x, y, xs, ys] = genNodesPositions(rAP, M, U)
% Inputs
%  rAP - AP coverage radius
%  M - number of APs
%  U - number of UTs
% Outputs
% [x, y] - APs positions
% [xs, ys] - UTs positions

% Generate APs positions
[x, y] = genUnifUserDist(M, rAP, 0, 0);

% Generate UTs positions within AP radius (500 m)
[xs, ys] = genUnifUserDist(U, rAP, 0, 0);
end

%% Generate uniform user distribution positions within a circle
function [x, y] = genUnifUserDist(M, radius, xCenter, yCenter)
r = radius * rand(1, M);
angles = 2 * pi + (0 - 2 * pi) .* rand(2, M);
x = xCenter + r .* cos(angles(1, :));
y = yCenter + r .* sin(angles(1, :));
end

function ApAng = genAngle()
ApAng = (rand - 0.5) * 180;
end

%% compute distance between two points
function d = mdist(A, B)
d = sqrt((A(1) - B(1))^2 + (A(2) - B(2))^2);
end

function [Ar, At] = genRadarChannel(c0, fc, P, Q, AngT, AngR)
elemT = (0:P-1)';
elemR = (0:Q-1)';

wavelen = c0 / fc;
d = wavelen / 2;
k = 2 * pi / wavelen;

At = exp(1j * elemT * k * d * sind(AngT)) / sqrt(P);
Ar = exp(1j * elemR * k * d * sind(AngR)) / sqrt(Q);
end

% Generate Truncated Laplace random variable
function y = lapRand(m, n, mu)
a = exp(-pi / 2); % angles lower than -pi/2 and higher than pi/2 are not considered
u1 = rand(m, n) - 0.5;
u2 = (1 - a) * (1 - 2 * abs(u1)) + a;

y = mu - 2 * sign(u1) .* log(u2);
y = mod(y + pi, 2 * pi) - pi; % constrain the angles to the range [-pi/2, pi/2]
end

