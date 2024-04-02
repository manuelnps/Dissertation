%function main
%% parameters
% User terminal parameters
U = 4;                              % Number of UEs
Ntx = 16;                           % Number of antennas UE

% Base station parameters
M = 2;                              % Number of BSs
Nrx = 16;                           % Number of antennas BS
NrxRF = U;                          % Number of RF chains BS

% Channel parameters
Ncl = 5;                            % Number of clusters
Nray = 3;                           % Number of rays
AngleSpread = 8;                    % Angle spread

% OFDM parameters
K = 64;                             % Number of subcarriers
MODULATION_ORDER = 4;               % QPSK

% Simulation parameters
Nchannels = 100;                    % Number of channel realizations
Nsym = 1;                           % Number of noise realizations
EbN0 = -20:10:20;                   % Energy per bit to noise power spectral density ratio
nvar = 10.^(-EbN0./10)/(log2(MODULATION_ORDER));   % Noise variance

%% Allocate memory and construct objects
ber = nan(Nchannels, Nsym, length(EbN0));
channel = CWideband_mmWave_Channel(Nrx, Ntx, K, K/4, Ncl, Nray, 10, 10, AngleSpread); % Construct channel object
DFT_S_OFDM = CDFT_S_OFDM(K, M, 1, MODULATION_ORDER);                             % Construct DFT-S-OFDM object

%% Open new figure to plot the results
figure('Position', [100, 100, 1*560, 1*420]);

%% Start timer 
tic
    
%% For each channel realization
for ch = 1:Nchannels
    
    % Generate channel for all BS and UE pairs
    H0 = generateChannels(channel, U, M);
    
    %% For each noise realization
    for n = 1:Nsym
        % Generate noise with unit variance
        noise = sqrt(1/2)*complex(randn(Nrx, K), randn(Nrx, K));
        
        %% Sweep over the EbN0 range
        for p = 1:length(EbN0)
            %% BS
            b = DFT_S_OFDM.genRandBits();                                       % Generate random bits
            s = DFT_S_OFDM.mod(b);                                              % DFT_S_OFDM modulation
            W = exp(1j*2*pi*rand(Ntx, M));                                      % Analog precoder
            
            %% Equivalent channel, combining channel and precoder
            H = getEquivalentChannel(H0, W);
            
            %% Channel
            y = zeros(Nrx, K, U);                                               % Allocate memory for received signal
            % Pass transmitted signal through the channel and add noise
            for u = 1:U
                for k = 1:K
                    y(:, k, u) = H(:, :, k, u)*s(k, :).' + sqrt(nvar(p))*noise(:,k);
                end
            end
            
            %% UE
            Ga = exp(1j*2*pi*rand(Nrx, NrxRF, U));                              % Random analog equalizer part           
            Gd = computeDigitalEqualizer(Ga, H, nvar(p));                       % MMSE digital equalizer part
            
            % Equalize received signal at each UE
            ce = zeros(M, K, U);
            for u = 1:U
                for k = 1:K
                    ce(:, k, u) = Gd(:, :, k, u)'*Ga(:, :, u)'*y(:, k, u);
                end
            end
            
            % Average signals from all UEs at BS
            ce = mean(ce, 3);
            
            % DFT-S-OFDM demodulation
            br = DFT_S_OFDM.demod(ce.');
            
            %% Compute BER
            [~, ber(ch, n, p)] = biterr(b, br);
        end
    end
    sim_time = toc; % Read timer
    
    %% Plot results (average BER over channel and noise)
    semilogy(EbN0, squeeze(mean(ber, [1 2])));
    axis([min(EbN0) max(EbN0) 1e-4 0.25])
    xlabel('EbN0 (dB)')
    ylabel('BER')
    pause(1)                                                                    % Wait for the plotting function to finish
    
    %% Display information about simulation
    if(mod(ch, 10) == 0)
        avg_sim_time = sim_time/ch;
        Nchannel_left = Nchannels - ch;
        fprintf('Simulation time left (estimate): %0.1f min. \n', avg_sim_time*Nchannel_left/60)
    end
end
%end


%% Generate equivalent channel, including precoder, for all BS and UE pairs
function H0 = generateChannels(channel, U, M)
% Channel parameters
np = 4.1;   % Path-Loss exponents
dp_s = 7.6; % Standard deviation of shadowing factor
f = 28;     % GHz
rUE = 50;   % UE coverage radius

% Allocate memory for variables
H0 = zeros(channel.Nr, channel.Nt, channel.Nc, U, M);
Arx0 = zeros(channel.Nr, channel.Ncl*channel.Nray, U, M);
pl = zeros(M, U);
muTx = zeros(M, channel.Ncl, U);

% Generate BS and UE positions
[xbs, ybs, xue, yue] = genNodesPositions(rUE, M, U);

% Average path loss for a scenario with a UE at cell border being served by a BS at cell center
averageShadowing = exp((log(10)/10*dp_s)^2/2);  % Theoretical value of average shadowing
pathLossWithoutShadowing = pLoss(f,rUE,np,0);   % Path loss not considering shadowing
C = pathLossWithoutShadowing*averageShadowing;  % Normalization constant

% For each BS
for m = 1:M
    % For each UE
    for u = 1:U
        % Generate channel without path loss
        [H0(:, :, :, u, m), ~, Arx0(:, :, u, m), ~, ~, muTx(m,:,u)] = channel.genRandomInstance;
        
        % Distance (D) between UEs and BSs:
        D = mdist([xue(u), yue(u)], [xbs(m), ybs(m)]);
        
        % Path Loss
        pl(m,u) = pLoss(f, D, np, dp_s)/C;
        
        H0(:, :, :, u, m) = sqrt(pl(m, u))*H0(:, :, :, u, m);
    end
end
end

%% Get equivalent channel combining channel and precoder
function H = getEquivalentChannel(H0, W)
G = H0(:,1,:,:,:);
for u = 1:size(H0, 4)
    for m = 1:size(H0, 5)
        for k = 1:size(H0, 3)
            G(:, :, k, u, m) = H0(:, :, k, u, m)*W(:, m);
        end
    end
end
H1 = G;

% Permute TX antenna with BS dim, as nr of TX antennas is one
H = permute(H1, [1 4 3 5 2]);
end

%% Compute digital part of equalizer
function Gd = computeDigitalEqualizer(Ga, H, nvar)
% Get parameters
Nrx = size(H, 1);                                   % Nr of UEs
U = size(H, 4);                                     % Nr of UEs
K = size(H, 3);                                     % Nr of subcarriers
M = size(H, 5);                                     % Nr of BSs
Nrf = size(Ga, 2);                                  % Nr of RF chains

% Allocate memory
Gd = zeros(Nrf, M, K, U);

% For each BS
for m = 1:M
    % Initialize auxiliary variable to be used to compute normalizing constant
    S0 = zeros(U, 1);
    % For each subcarrier
    for k = 1:K
        % Compute received signal correlation
        R = H(:, :, k, :, m)*H(:, :, k, :, m)' + nvar*eye(Nrx);
        % Compute optimum digital part
        % Use the pseudo-inverse instead of inverse since the matrix R may
        % be ill-conditioned as some UE-BS pairs may have very low path loss
        Gd(:,:,k,m) = inv(Ga(:,:,m)'*R*Ga(:,:,m))*(Ga(:,:,m)'*H(:,:,k,:,m));
        % Compute auxiliary random variable to normalize digital part
        S0 = S0 + diag(Gd(:,:,k,m)'*Ga(:,:,m)'*H(:,:,k,:,m));
    end
    
    % For each subcarrier
    for k = 1:K
        % Normalize digital part such that \sum_k diag(Gd,k'*Ga'*Hk) = I
        Gd(:,:,k,m) = Gd(:,:,k,m)*diag(K./S0.'');
    end
end
end

%% Path loss
function PL= pLoss(f,d,E_NLOS,f_NLOS)
waveLen = 0.3/f;

G_TX=15;      % 15 dBi
G_Rx=24.5;    % 24.5 dBi
Gain=G_TX+G_Rx;
Gain_L=10^(Gain/10);

d0=1;
Beta0=10*log10((4*pi*d0)/waveLen).^2;

%E_NLOS=4.1;  % NLOS Path-loss exponents in dB
%f_NLOS=7.6;  % NLOS Standard deviation of shadowing factor in dB

A_NLOS=f_NLOS*randn();     % A_NLOS= 10.^(f_NLOS*randn()/10);  % Linear

Beta= Beta0 + 10*E_NLOS*log10(d/d0) + A_NLOS;
Beta_L= 10^(Beta/10);

% Ga_Beta= Gain - Beta;
% Ga_Beta_L1=10^(Ga_Beta/10);

PL= Gain_L/Beta_L;
end

%% Generate nodes positions
function [x, y, xs, ys] = genNodesPositions(rUE, M, U)
% Inputs
% rUE - UE coverage radius
% M - Number of BSs
% U - Number of UEs
% Outputs
% [x, y] - BSs positions
% [xs, ys] - UEs positions

%%%%%%%%%%%%%% Code %%%%%%%%%%%%%
% Generate BSs positions
[x, y] = genUnifUserDist(M, rUE, 0, 0);

% Generate UEs positions within BS radius (50 m)
[xs, ys] = genUnifUserDist(U, rUE, 0, 0);
end

%% Generate uniform user distribution positions within a circle
function [x, y] = genUnifUserDist(M, radius, xCenter, yCenter)
r = radius*rand(1, M);
angles = 2*pi + (0-2*pi).*rand(2, M);
x = xCenter + r.*cos(angles(1,:));
y = yCenter + r.*sin(angles(1,:));
end

%% Compute distance between two points
function d = mdist(A, B)
d = sqrt((A(1) - B(1)).^2 + (A(2) - B(2)).^2);
end
