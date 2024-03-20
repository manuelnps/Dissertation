% wideband mmWave channel model class
classdef CWideband_mmWave_Channel < handle
    properties
        fc {mustBePositive, mustBeReal} = 72;               % carrier frequency (in GHz)
        AngleSpread {mustBeNonnegative, mustBeReal} = 10;   % Angle spread in degrees
        dp_phi {mustBePositive, mustBeReal}                 % Laplace dist. scale parameter for the selected angle spread
        Ncl {mustBePositive, mustBeInteger}                 % number of clusters
        Nray {mustBePositive, mustBeInteger}                % number of rays
        Nt {mustBePositive, mustBeInteger}                  % number of transmit antennas
        Nr {mustBePositive, mustBeInteger}                  % number of receive antennas
        Nc {mustBePositive, mustBeInteger}                  % number of subcarriers
        Ncp {mustBePositive, mustBeInteger}                 % cyclic prefix length
        decCl {mustBeReal} = 10;                            % power decay clusters (in dB)
        decRa {mustBeReal} = 10;                            % power decay rays (in dB)
        
        % pulse shapping function parameters
        beta {mustBeReal} = 1;                              % roll-off factor
        samples {mustBePositive, mustBeInteger} = 1;        % nr samples per symbol
        pulseCoef;                                          % pulse coefficients
    end
    methods
        % constructor
        function obj = CWideband_mmWave_Channel(Nr, Nt, Nc, Ncp, Ncl, Nray, decCl, decRa, AngleSpread)
            obj.dp_phi = getlapRandParameter(obj.AngleSpread*pi/180);
            obj.Ncl = Ncl;          obj.Nray = Nray;
            obj.Nt = Nt;            obj.Nr = Nr;
            obj.Nc = Nc;            obj.Ncp = Ncp;
            obj.decCl = decCl;      obj.decRa = decRa;
            obj.AngleSpread = AngleSpread;
            obj.pulseCoef = rcosdesign(obj.beta,2*obj.Ncp,obj.samples);
        end
        
        % function to generate random channel matrix
        function [H, h, Arx, Atx, avgAoA, avgAoD] = genRandomInstance(obj)
            waveLen = 0.3/obj.fc;
            d = waveLen/2;
            k = 2*pi/waveLen;
            
            decTempCl = obj.Ncl*obj.Nray/log( 10^(obj.decCl/10) );
            t = 0:obj.Ncl*obj.Nray-1;
            kcl = exp(-t/decTempCl);
            kcl = kcl/sqrt(sum(kcl.^2));
            kcl = kcl*sqrt(obj.Ncl*obj.Nray );

            decTempRay = obj.Nray/log( 10^(obj.decRa/10) );
            t = 0:obj.Nray-1;
            kra = exp(-t/decTempRay);
            kra = kra/sqrt(sum(kra.^2));
            kra = kra*sqrt( obj.Nray );
            kra = repmat(kra,1,obj.Ncl);

            gama = sqrt((obj.Nt*obj.Nr)/(obj.Ncl*obj.Nray));

            elemT = (0:obj.Nt-1).';
            elemR = (0:obj.Nr-1).';

            passoMedio = floor( obj.Ncp /(obj.Ncl*obj.Nray)); tp=1; nra=1;

            avgAoD = zeros(obj.Ncl, 1);
            avgAoA = zeros(obj.Ncl, 1);
            
            A = zeros(obj.Nr*obj.Nt, obj.Ncl*obj.Nray);
            yp = zeros(obj.Ncl*obj.Nray, obj.Ncp);
            
            Atx = zeros(obj.Nt, obj.Ncl*obj.Nray); 
            Arx = zeros(obj.Nr, obj.Ncl*obj.Nray);
            for cl=1:obj.Ncl
                avgAoD(cl) = 2*pi*rand() - pi;
                avgAoA(cl) = 2*pi*rand() - pi;
                for ray=1:obj.Nray
                    % Generate Rays with Truncated Laplace Distribution
                    fasesT = lapRand(1,1, avgAoD(cl), obj.dp_phi);
                    At = exp(1j*elemT*k*d*sin(fasesT))/sqrt(obj.Nt);
                    Atx(:, (cl - 1)*obj.Nray + ray) = At;

                    % Generate Rays with Truncated Laplace Distribution
                    fasesR = lapRand(1,1, avgAoA(cl), obj.dp_phi);
                    Ar = exp(1j*elemR*k*d*sin(fasesR))/sqrt(obj.Nr);
                    Arx(:, (cl - 1)*obj.Nray + ray) = Ar;

                    alfa = (randn() + 1j*randn())/sqrt(2);

                    yp(nra, :) = kcl(nra)*kra(nra)*gama*alfa*p(tp, obj.pulseCoef, obj.Ncp, obj.samples);

                    if passoMedio==1
                        tp = nra*passoMedio;
                    else
                        tp = nra*passoMedio + randi([-floor(passoMedio/2)+1 floor(passoMedio/2)-1]);
                    end
                    nra = nra+1;
                end
            end
            
            % The following is equivalent to
            % for m = 1:obj.Ncp
            %     h(:,:,m) = Arx*diag(yp(:, m))*Atx';
            % end
            for m = 1:obj.Ncl*obj.Nray
                tmp = Arx(:, m)*Atx(:, m)';
                A(:, m) = tmp(:);
            end
            h = reshape(A*yp, obj.Nr, obj.Nt, obj.Ncp);
            
            H = fft(h, obj.Nc, 3);
        end
    end
end

% window function
function y = p(t, h, Ncp, samples)
kmax = samples*Ncp+1;
y = h(kmax-t+1:kmax-t+Ncp);
y = y/sum(y);
end

% get Laplace distribution b parameter
function x = getlapRandParameter(sigma)
xlow = 0; xhigh = pi;
for n = 1:100
    x = (xlow + xhigh)/2;
    if(2*x^2 - pi*(2*x + pi)/(exp(pi/x) - 1) < sigma^2)
        xlow = x;
    else
        xhigh = x;
    end
end
end

% Generate Truncated Laplace random variable
function y  = lapRand(m, n, mu, b)
a = exp(-pi/b);   % angles lower than -pi and higher than pi are not considered
u1 = rand(m, n) - 0.5;
u2 = (1 - a)*(1 - 2*abs(u1)) + a;

y = mu - b * sign(u1).* log(u2);
y = mod(y + pi, 2*pi) - pi; % contraint the angles to the range [-pi, pi]
end







    