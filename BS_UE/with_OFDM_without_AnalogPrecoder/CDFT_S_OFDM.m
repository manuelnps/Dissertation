%
% The Discrete Fourier Transform Spread (DFT-S) is a technique used in wireless communication systems, 
% particularly in Orthogonal Frequency Division Multiplexing (OFDM) systems, to achieve certain desirable properties such as 
% spreading the spectrum of the transmitted signal and providing robustness against frequency-selective fading channels.
% 
% In DFT-S, spreading is achieved by multiplying the symbols to be transmitted by the Discrete Fourier Transform (DFT) matrix. 
% This spreads the symbols across multiple subcarriers, effectively increasing the bandwidth occupied by the signal. 
% This spreading operation is typically followed by an Inverse Fast Fourier Transform (IFFT) to convert the signal into the time 
% domain before transmission.
% 
% The spreading process serves multiple purposes:
% 
% Spectrum Spreading: By spreading the symbols across multiple subcarriers, the signal occupies a wider bandwidth. 
% This can help combat narrowband interference and improve robustness against frequency-selective fading channels.
% 
% Orthogonality: The DFT matrix used for spreading ensures that the subcarriers are orthogonal to each other. 
% This orthogonality property is crucial for mitigating inter-symbol interference (ISI) and simplifying equalization at the receiver.
% 
% Frequency Diversity: By spreading the signal across multiple subcarriers, DFT-S-OFDM systems exploit frequency diversity. 
% This means that even if certain subcarriers experience deep fades due to multipath propagation, other subcarriers may still 
% carry useful information, improving the overall reliability of the system.

% DFT-S-OFDM class
classdef CDFT_S_OFDM
    properties
        Nc {mustBePositive, mustBeInteger};                                 % number of subcarriers
        Ns {mustBePositive, mustBeInteger} = 1;                             % Number of streams (number of independent data streams)
        Ncp {mustBePositive, mustBeInteger};                                % CP length
        SF {mustBePositive, mustBeInteger} = 4;                             % Spreading factor
        MODULATIONORDER {mustBePositive, mustBeInteger} = 4;                % modulation order
        pr {mustBeReal};                                                    % random permutation to map subcarriers
    end
    methods
        %% classs constructor
        function obj = CDFT_S_OFDM(Nc, Ns, SF, MODULATIONORDER)
            obj.Nc = Nc;        % set nr subcarriers
            obj.Ns = Ns;        % set nr streams
            obj.SF = SF;        % set spreading factor
            obj.Ncp = 16;
            obj.MODULATIONORDER = MODULATIONORDER;  % set MODULATIONORDER
            obj.pr = randperm(Nc);  % generate random permutation
        end
        
        %% generate random bits
        function b = genRandBits(obj)
            b = randi([0, 1], obj.Nc*log2(obj.MODULATIONORDER), obj.Ns);
        end
        
        %% modulation
        % generate random bits; modulate bits using M-QAM; apply DFT-spreading
       function  s = mod(obj, b)
            % M-QAM modulation
            s = qammod(b, obj.MODULATIONORDER, 'gray', 'InputType', 'bit', 'UnitAveragePower', true);
            
            
        end
        
        %% demodulation
        function b = demod(obj, s)
           
            % M-QAM demodulation
            b = qamdemod(s, obj.MODULATIONORDER, 'gray', 'OutputType', 'bit', 'UnitAveragePower', true);
        end
        
        %% DFT-S-OFDM spreading
        function y = DFTspread(obj, x)
            y = reshape(fft(reshape(x, obj.Nc/obj.SF, obj.SF, obj.Ns), [], 1)/sqrt(obj.Nc/obj.SF), obj.Nc, obj.Ns);
            y = y(obj.pr, :);
        end
        
        %% DFT-S-OFDM unspreading
        function y = DFTunspread(obj, x)
            x(obj.pr, :) = x;
            y = reshape(ifft(reshape(x, obj.Nc/obj.SF, obj.SF, obj.Ns), [], 1)*sqrt(obj.Nc/obj.SF), obj.Nc, obj.Ns);
        end
    end
end