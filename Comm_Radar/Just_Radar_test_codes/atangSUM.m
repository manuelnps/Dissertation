clc;
clear;

% Parameters
M = 4;
Ntargs = 3; % Number of targets
c0 = 3e8;
fc = 28e9; % Carrier frequency
lambda = c0 / fc;
NAp = 16;
elemT = (0:NAp-1)';
elemR = (0:NAp-1)';
wavelen = c0 / fc;
d = wavelen / 2;
k = 2 * pi / wavelen;

% Generate random angles
ApAng = (rand - 0.5) * 180;
AngT = (rand(1, M) - 0.5) * 180;
At = exp(1j * elemT * k * d * sind(AngT)) / sqrt(NAp);
AtH = At';

AngR = (rand(1, Ntargs) - 0.5) * 180;
Ar = exp(1j * elemR * k * d * sind(AngR)) / sqrt(NAp);

% Received signal
x = zeros(NAp, M, Ntargs);
for n = 1:Ntargs
    for m = 1:M
        x(:,m,n) = Ar(:,n) * (AtH(m,:) * At(:,m)); 
    end 
end

% Summing received signals
RxSignal = sum(sum(x, 3), 2);

% Angle estimation using periodogram
angles = -90:0.1:90; % Angle grid
P = zeros(size(angles));

for idx = 1:length(angles)
    steering_vector = exp(1j * elemR * k * d * sind(angles(idx))) / sqrt(NAp);
    P(idx) = abs(steering_vector' * RxSignal)^2;
end

% Normalize periodogram
P = P / max(P);

% Plot periodogram
figure;
plot(angles, P);
xlabel('Angle (degrees)');
ylabel('Normalized Power');
title('Periodogram for Angle Estimation');
grid on;

% Display estimated angle
[~, maxIdx] = max(P);
estimated_angle = angles(maxIdx);
disp(['Estimated Angle: ', num2str(estimated_angle), ' degrees']);