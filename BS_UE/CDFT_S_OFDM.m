% DFT-S-OFDM class
classdef CDFT_S_OFDM
    properties
        Nc {mustBePositive, mustBeInteger};                                 % number of subcarriers
        Ns {mustBePositive, mustBeInteger} = 1;                             % Number of streams
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
            obj.MODULATIONORDER = MODULATIONORDER;  % set MODULATIONORDER
            obj.pr = randperm(Nc);  % generate random permutation
        end
        
        %% generate random bits
        function b = genRandBits(obj)
            b = randi([0, 1], obj.Nc*log2(obj.MODULATIONORDER), obj.Ns);
        end
        
        %% modulation
        % generate random bits; modulate bits using M-QAM; apply DFT-spreading
        function [c, s] = mod(obj, b)
            % M-QAM modulation
            s = qammod(b, obj.MODULATIONORDER, 'gray', 'InputType', 'bit', 'UnitAveragePower', true);
            
            % DFT spreading
            c = DFTspread(obj, s);
        end
        
        %% demodulation
        function [b, s] = demod(obj, c)
            % DFT unspreading
            s = DFTunspread(obj, c);
            
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