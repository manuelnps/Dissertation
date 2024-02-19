function main
%% parameters
% user terminal parameters
U = 2;                              % nr Users
Ntx = 16;                           % nr antennas UT

% AP parameters
M = 4;                              % nr APs
Nrx = 16;                           % nr antennas AP
NrxRF = U;                          % nr RF chains AP

% channel parameters
Ncl = 5;                            % nr clusters
Nray = 3;                           % nr rays
AngleSpread = 8;                    % angle spread

% OFDM parameters
K = 64;                            % nr of subcarriers
MODULATION_ORDER = 4;               % QPSK



% simulation parameters
Nchannels = 100;                    % number of channel realizations
Nsym = 1;                          % number of noise realizations
EbN0 = -20:10:20;                   % energy per bit to noise power spectral density ratio
nvar = 10.^(-EbN0./10)/(log2(MODULATION_ORDER));   % noise variance

%% allocate memory and construct objects
ber = nan(Nchannels, Nsym, length(EbN0));
channel = CWideband_mmWave_Channel(Nrx, Ntx, K, K/4, Ncl, Nray, 10, 10, AngleSpread);%% construct channel object
DFT_S_OFDM = CDFT_S_OFDM(K, U, 1, MODULATION_ORDER);                            %% construct DFT-S-OFDM object

%% open new figure to plot the results
figure('Position', [100, 100, 1*560, 1*420]);

%% start timer 
tic
    
%% for each channel realization
for ch = 1:Nchannels
    
    % generate channel for all UT and AP pairs
    H0 = generateChannels(channel, M, U);
    
    %% for each noise realization
    for n = 1:Nsym
        % generate noise with unit variance
        noise = sqrt(1/2)*complex(randn(Nrx, K), randn(Nrx, K));
        
        %% sweep over the EbN0 range
        for p = 1:length(EbN0)
            %% User terminals
            b = DFT_S_OFDM.genRandBits();                                       % Generate random bits
            s = DFT_S_OFDM.mod(b);                                              % DFT_S_OFDM modulation
            W = exp(1j*2*pi*rand(Ntx, U));                                      % analog precoder
            
            %% Equivalent channel, combining channel and precoder
            H = getEquivalentChannel(H0, W);
            
            %% Channel
            y = zeros(Nrx, K, M);                                               % alocate memory for received signal
            % pass transmitted signal through the channel and add noise
            for m = 1:M
                for k = 1:K
                    y(:, k, m) = H(:, :, k, m)*s(k, :).' + sqrt(nvar(p))*noise(:,k);
                end
            end
            
            %% Base station
            Ga = exp(1j*2*pi*rand(Nrx, NrxRF, M));                              % random analog equalizer part           
            Gd = computeDigitalEqualizer(Ga, H, nvar(p));                       % MMSE digital equalizer part
            
            % equalize received signal at each AP
            ce = zeros(U, K, M);
            for m = 1:M
                for k = 1:K
                    ce(:, k) = Gd(:, :, k, m)'*Ga(:, :, m)'*y(:, k, m);
                end
            end
            
            % average signals from all APs at CU
            ce = mean(ce, 3);
            
            % DFT-S-OFDM demodulation
            br = DFT_S_OFDM.demod(ce.');
            
            %% compute ber
            [~, ber(ch, n, p)] = biterr(b, br);
        end
    end
    sim_time = toc; % read timer
    
    %% plot results (average ber over channel and noise)
    semilogy(EbN0, squeeze(nanmean(ber, [1 2])));
    axis([min(EbN0) max(EbN0) 1e-4 0.25])
    xlabel('EbN0 (dB)')
    ylabel('BER')
    pause(1)                                                                    % wait for the plotting function to finish
    
    %% display information about simulation
    if(mod(ch, 10) == 0)
        avg_sim_time = sim_time/ch;
        Nchannel_left = Nchannels - ch;
        fprintf('Simulation time left (estimate): %0.1f min. \n', avg_sim_time*Nchannel_left/60)
    end
end
end


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

%% get equivalent channel combining channel and precoder
function H = getEquivalentChannel(H0, W)
G = H0(:,1,:,:, :);
for m = 1:size(H0, 4)
    for u = 1:size(H0, 5)
        for k = 1:size(H0, 3)
            G(:, :, k, m, u) = H0(:, :, k, m, u)*W(:, u);
        end
    end
end
H1 = G;

% permute tx antenna with user dim, as nr of tx antennas is one
H = permute(H1, [1 5 3 4 2]);
end
%% compute digital part of equalizer
function Gd = computeDigitalEqualizer(Ga, H, nvar)
%% get parameters
Nrx = size(H, 1);                                   % nr of users
U = size(H, 2);                                     % nr of users
K = size(H, 3);                                     % nr of subcarriers
M = size(H, 4);                                     % nr of APs
Nrf = size(Ga, 2);                                  % nr of RF chains

%% allocate memory
Gd = zeros(Nrf, U, K, M);

% for each AP
for m = 1:M
    % initialize auxiliary variable to be used to compute normalizing constant
    S0 = zeros(U, 1);
    % for each subcarrier
    for k = 1:K
        % compute received signal correlation
        R = H(:, :, k, m)*H(:, :, k, m)' + nvar*eye(Nrx);
        % compute optimum digital part
        % is used the pseudo inverse instead of inverse since the matrix R may
        % be ill conditioned as some UT-AP pairs may have very low path loss
        Gd(:,:,k,m) = inv(Ga(:,:,m)'*R*Ga(:,:,m))*(Ga(:,:,m)'*H(:,:,k,m));
        % compute auxiliary random variable to normalize digital part
        S0 = S0 + diag(Gd(:,:,k,m)'*Ga(:,:,m)'*H(:,:,k,m));
    end
    
    % for each subcarrier
    for k = 1:K
        % normalize digital part such that \sum_k diag(Gd,k'*Ga'*Hk) = I
        Gd(:,:,k,m) = Gd(:,:,k,m)*diag(K./S0.'');
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

