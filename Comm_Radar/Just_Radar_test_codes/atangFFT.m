clc;
clear;

% Parameters
M = 16;                  % Number of transmit antennas
Ntargs = 1;             % Number of targets
c0 = 3e8;               % Speed of light in m/s
fc = 28e9;              % Carrier frequency in Hz
lambda = c0 / fc;       % Wavelength
NAp = 16;               % Number of array elements
elemT = (0:NAp-1)';     % Transmit array element indices
elemR = (0:NAp-1)';     % Receive array element indices
d = lambda / 2;         % Element spacing
k = 2 * pi / lambda;    % Wave number

% Generate random angles for targets and antenna arrays
ApAng = (rand - 0.5) * 180;
AngT = (rand(1, M) - 0.5) * 180;
At = exp(1j * elemT * k * d * sind(AngT)) / sqrt(NAp); % Transmit steering vectors
AtH = At'; % Hermitian transpose of At

%AngR = (rand(1, Ntargs) - 0.5) * 180;
AngR = 30;
Ar = exp(1j * elemR * k * d * sind(AngR)) / sqrt(NAp); % Receive steering vectors

% Generate received signal matrix
x = zeros(NAp, M, Ntargs);
for n = 1:Ntargs
    for m = 1:M
        x(:,m,n) = Ar(:,n) * (AtH(m,:) * At(:,m))*randi(1,1); 
    end 
end

% Summing received signals across targets and antennas
RxSignal = sum(sum(x, 3), 2);

% Compute periodogram using one FFT
N = length(RxSignal); % Length of RxSignal
angles = 0:180; % Angle grid for estimation
P = zeros(size(angles));

% FFT of the received signal
S = fft(RxSignal, N);

% Precompute steering vectors for all angles
steering_vectors = zeros(NAp, length(angles));
for idx = 1:length(angles)
    angle = angles(idx);
    steering_vectors(:, idx) = exp(1j * elemR * k * d * sind(angle-91)) / sqrt(NAp);
end

% FFT of all steering vectors
W = fft(steering_vectors, N); %%nao e igual

% Compute periodogram using cross-power spectral density approach
for idx = 1:length(angles)
    P(idx) =( 1/ N)*abs(W(:, idx)' * S)^2;
end

% Normalize periodogram
P = P / max(P);
xaxis = -90:90;
% Plot periodogram
figure;
plot(xaxis, P);
xlabel('Angle (degrees)');
ylabel('Normalized Power');
title('Periodogram for Angle Estimation (FFT)');
grid on;
