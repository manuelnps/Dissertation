function OFDM_Radar
clear all
clc
%% add all constants here
% Waveform parameters

fc=24e9;  %carrier frequency
B=93.1e6; %signal bandwidth
M=256; %Number of OFDM symbols
N=32;  %Number of sub-carriers

 % Sub-carrier separation
Tsymbol = N/B ; % Useful symbol duration
Tcp=(1/4)*Tsymbol;
To=Tsymbol+Tcp;

spacing=1/Tsymbol;

% fc=1.2932e13; %carrier frequency
% B=93.1e6; %signal bandwidth
% M=16; %Number of OFDM symbols
% N=8;  %Number of sub-carriers
% spacing = 9.3168e+07; %15e3; % Sub-carrier separation = 15 kHz
% Tsymbol = 1/spacing; % Useful symbol duration = 
% Tcp = 5.21e-6  ; % Cyclic prefix duration = 5.21 us
% To= Tsymbol + Tcp; %OFDM symbol duration = 71.86 us
% 
% Ts=1/(spacing*N); %sampling time
% fs=1/Ts;



%Target Parameters
c0=3e8;
%Resolution and Unambiguous Range/Velocity
resolution_R=c0/(2*spacing*N); %Resolution Range
resolution_v=c0/(2*fc*Tsymbol*M);%Resolution Velocity
%R= m , v= m/s
d=0.5; %Distance m
v=0.3; %velocity m/s

delay=2*d/c0; %delay
fd=2*(v/c0)*fc; %Doppler Shift

%Periodogram
Nper=8*N;
Mper=8*M;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%Transfer function F
F=zeros(N,M);

for k=1:N
   for l=1:M
        F(k,l)=exp(j*2*pi*(l-1)*To*fd)*exp(-j*2*pi*(k-1)*delay*spacing);
   end    
end


F_final= F.'; %Matrix Number of OFDM symbols by Number of sub-carriers


% Two-dimension Periodogram
Per=zeros(Nper,Mper);
Aux=zeros(N,Mper);


Aux=fft(F.',Mper).';
%Aux=fftshift(Aux.').';
Per=ifft(Aux,Nper);
Per=fftshift(Per);
Per=(1/(N*M)).*abs(Per).^2;

Per_final=Per.';


%Plot
y=[(-Nper/2):((Nper/2)-1)].*(c0/(2*spacing*Nper));
x=[(-Mper/2):((Mper/2)-1)].*(c0/(2*fc*To*Mper));

figure(1);
contour(x,y,Per)
title('Periodogram')
xlabel('Relative speed m/s')
ylabel('Distance m')


figure(2);
imagesc(x,y,Per)
title('Periodogram')
xlabel('Relative speed m/s')
ylabel('Distance m')
colorbar


end