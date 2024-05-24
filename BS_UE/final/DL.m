% fazer com mais parametros diferentes n aps, n UE

load('4U16AP.mat');
% Plot the graph
figure('Position', [100, 100, 1*560, 1*420]);
title('UE - AP');
semilogy(EbN0, squeeze(mean(ber, [1 2])), 'LineWidth', 2,'DisplayName','UE - AP');
axis([min(EbN0) max(EbN0) 1e-4 0.65]);
xlabel('EbN0 (dB)');
ylabel('BER');
grid on;
hold on;
clc;
clear;


%% For comparing with UE - AP data
clc;
clear;
%function main
%% parameters
% AP parameters
M = 16;                              % n APs
NAp = 16;                           % nr antennas Acess Point

% UE  parameters
U = 4;                              % nr UEs
NUe = 1;                           % nr antennas User Equipment
NrxRF = NUe;                        % nr RF chains is equal to the number of antennas

% channel parameters


Ncl = 5;                            % nr clusters
Nray = 3;                           % nr rays
AngleSpread = 8;                    % angle spread

% OFDM parameters
K = 64;                            % nr of subcarriers
MODULATION_ORDER = 4;               % QPSK



% simulation parameters
Nchannels = 1000;                    % number of channel realizations
Nsym = 1;                          % number of noise realizations
EbN0 = -28:4:4;                   % energy per bit to noise power spectral density ratio 
%EbN0=100;
nvar = 10.^(-EbN0./10)/(log2(MODULATION_ORDER));   % noise variance

%% allocate memory and construct objects
ber = nan(Nchannels, Nsym, length(EbN0));
channel = CWideband_mmWave_Channel(NUe, NAp, K, K/4, Ncl, Nray, 10, 10, AngleSpread);%% construct channel object (DL)
DFT_S_OFDM = CDFT_S_OFDM(K, U, 1, MODULATION_ORDER);                            %% construct DFT-S-OFDM object //sc-fdma (DL)

%% open new figure to plot the results

%% start timer 
tic
    
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
            
            Wd = computeDigitalPrecoder(H, nvar(p),M);

            y1 = zeros(K,U);
            y = zeros(K,U); 
            
            for k = 1:K
                for u = 1:U
                    for m = 1:M
                        y1(k,u) = y1(k,u)+H0(:, :, k, u, m) * Wd(:, :, k, m) * s(k, :).' ;
                        
                    end
                    y(k,u) = y1(k,u) + sqrt(nvar(p))*noise(:,k);
                end
            end
     
            br = DFT_S_OFDM.demod(y);
           
            [~, ber(ch, n, p)] = biterr(b, br);
            
        end
    end
    sim_time = toc; % read timer
    
    %% plot results (average ber over channel and noise)
    %semilogy(EbN0, squeeze(nanmean(ber, [1 2])));
    %semilogy(EbN0, squeeze(mean(ber, [1 2])));
    figure(1)
    semilogy(EbN0, squeeze(mean(ber, [1 2])), 'r', 'LineWidth', 2, 'DisplayName','UE - AP');
    axis([min(EbN0) max(EbN0) 1e-4 0.65]); % Adjust the y-axis limit
    xlabel('EbN0 (dB)')
    ylabel('BER')
    grid on;
    sgtitle('Comparison between UE - AP and AP - UE');
    pause(1)                                                                    % wait for the plotting function to finish
    
    %% display information about simulation
    if(mod(ch, 10) == 0)
        avg_sim_time = sim_time/ch;
        Nchannel_left = Nchannels - ch;
        fprintf('Simulation time left (estimate): %0.1f min. \n', avg_sim_time*Nchannel_left/60)
    end
end
%end
legend({'UE-AP','AP-UE'})

%% generate equivalent channel, including precoder, for all UT and AP pairs
function H0 = generateChannels(channel, M, U)
%% Channel parameters
np = 4.1;   % path-Loss exponents
dp_s = 7.6; % standard deviation of shadowing factor
f = 28;     % GHz
rAP = 500;  % AP coverage radius

%% allocate memory for variables
H0 = zeros(channel.Nr, channel.Nt, channel.Nc, M, U);
Arx0 = zeros(channel.Nr, channel.Ncl*channel.Nray, M, U);
pl = zeros(U, M);
muTx = zeros(U, channel.Ncl, M);

%% generate users and access points positions
[xa, ya, xu, yu] = genNodesPositions(rAP, M, U);

%% Average path loss for a scenario with a user at cell border being served by an AP at cell center
averageShadowing = exp((log(10)/10*dp_s)^2/2);  % theoretical value of average shadowing
pathLossWithoutShadowing = pLoss(f,rAP,np,0);   % path loss not considering shadowing
C = pathLossWithoutShadowing*averageShadowing;  % normalization constant

% for each AP
for m = 1:M
    % for each UT
    for u = 1:U
        % generate channel without path loss
        [H0(:, :, :, m, u), ~, Arx0(:, :, m, u), ~, ~, muTx(u,:,m)] = channel.genRandomInstance;
        
        % distance (D) between UTs and APs:
        D = mdist([xu(u), yu(u)], [xa(m), ya(m)]);
        
        % Path Loss
        pl(u,m) = pLoss(f, D, np, dp_s)/C;
        
        H0(:, :, :, m, u) = sqrt(pl(u, m))*H0(:, :, :, m, u);
    end
end
end


function  Wd  = computeDigitalPrecoder(H, nvar, M)
    % Get parameters
    NAp = size(H, 2);  % Number of receivers
    U = size(H, 1);    % Number of users
    K = size(H, 3);    % Number of subcarriers
    M = size(H, 4);    % Number of APs
    
    % Allocate memory for Gd
    Wd = zeros(NAp, U, K, M);
    
    S0 = zeros(U, 1);
    % Compute Gd without Ga
    for m = 1:M
        
        for k = 1:K
            % Compute received signal correlation
            R = H(:, :, k, m) * H(:, :, k, m)'+ nvar * eye(U);
            % Compute optimum digital part
            Wd(:, :, k, m) =  H(:, :, k, m)'*inv(R);
            powers = diag(Wd(:,:,k,m)'*Wd(:,:,k,m));
            Wd(:,:,k,m) = Wd(:,:,k,m)/sqrt(M*diag(powers));
        end
    end
end


%% Path loss
function PL= pLoss(f,d,E_NLOS,f_NLOS)
waveLen = 0.3/f;

G_TX=15;      %15 dBi
G_Rx=24.5;    %24.5 dBi
Gain=G_TX+G_Rx;
Gain_L=10^(Gain/10);

d0=1;
Beta0=10*log10((4*pi*d0)/waveLen).^2;

%E_NLOS=4.1;  %NLOS Path-loss exponents in dB
%f_NLOS=7.6;  %NLOS Standard deviation of shadowing factor in dB

A_NLOS=f_NLOS*randn();     %A_NLOS= 10.^(f_NLOS*randn()/10);  %linear

Beta= Beta0 + 10*E_NLOS*log10(d/d0) + A_NLOS;
Beta_L= 10^(Beta/10);

%Ga_Beta= Gain - Beta;
%Ga_Beta_L1=10^(Ga_Beta/10);

PL= Gain_L/Beta_L;
end

%% Generate nodes positions
function [x, y, xs, ys] = genNodesPositions(rAP, M, U)
% Inputs
%  rAP - AP coverage radius
%  M - number of APs
%  U - number of UTs
% Outputs
%  [x, y] - APs positions
%  [xs, ys] - UTs positions

%%%%%%%%%%%%%% Code %%%%%%%%%%%%%
% Generate APs positions
[x, y] = genUnifUserDist(M, rAP, 0, 0);

% Generate UTs positions within AP radius (500 m)
[xs, ys] = genUnifUserDist(U, rAP, 0, 0);
end

%% Generate uniform user distribution positions within a circle
function [x, y] = genUnifUserDist(M, radius, xCenter, yCenter)
r = radius*rand(1, M);
angles = 2*pi + (0-2*pi).*rand(2, M);
x = xCenter + r.*cos(angles(1,:));
y = yCenter + r.*sin(angles(1,:));
end

%% compute distance between two points
function d = mdist(A, B)
d = sqrt((A(1) - B(1)).^2 + (A(2) - B(2)).^2);
end

