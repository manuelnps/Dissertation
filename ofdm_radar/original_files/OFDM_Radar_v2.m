function OFDM_Radar
clear all
clc
%% add all constants here
% Waveform parameters

fc=24e9;  %carrier frequency
B=93e6; %signal bandwidth
M=256; %Number of OFDM symbols
N=1024;  %Number of sub-carriers


% Sub-carrier separation
Tsymbol = N/B ; % Useful symbol duration
Tcp=(1/4)*Tsymbol;
To=Tsymbol+Tcp;

spacing=1/Tsymbol;

Ts = 1/(spacing*N); %Sampling time
fs = 1/Ts; %frequency sampling
NCP = round(Tcp/Ts);

data_mod2=zeros(1,N);

%Target Parameters
c0=3e8;
%Resolution and Unambiguous Range/Velocity
resolution_R=c0/(2*spacing*N); %Resolution Range
resolution_v=c0/(2*fc*Tsymbol*M);%Resolution Velocity

%H>1 
%R= m , v= m/s
d=[0 100]; %Distance m
v=[0 30]; %velocity m/s
delay=[];
fd=[];

for h=1:length(d)
    delay(h)=2*d(h)/c0; %delay
    fd(h)=(2*(fc/c0))*v(h); %Doppler Shift
end


%Periodogram
Nper=8*N;
Mper=8*M;

%Transfer Function
for h=1:length(delay)
Ftx=zeros(1,1);
s_ifft_comp=zeros(1,N);
F=zeros(N,M);
sum=0;
Frx_final=zeros(N,M);
Frx=zeros(N,M);
Frx_new=zeros(N,M);
%%%%%%%%%%%%%%%%%%%%%%%%%
%%
   
  
    for l=1:M

    %%%%Transmitter

        %Data
        nbits=N*log2(2); %BPSK
        data=round(rand(1,nbits)); %data generation
        data_mod = -data*2+1; %data modulation
        Ftx=[Ftx data_mod];


        %OFDM frame + CP

       
       %Apply delay H>1
       
       %data_mod2=data_mod2+data_mod.*exp(-2*j*pi*delay(h)*spacing*(0:(N-1)));
       data_mod2=data_mod.*exp(-2*j*pi*delay(h)*spacing*(0:(N-1)));
       % data_mod=data_mod.*exp(-2*1i*pi*delay(h)*spacing*(0:(N-1)));

        % Inverse FOURIER Transform
        s_ifft_comp=ifft(data_mod2);

        % Add CP

        OFDM_frame_cp2 = [s_ifft_comp(end-NCP+1:end) s_ifft_comp];
        %OFDM_frame_cp2 =[OFDM_frame_delay(end-NCP+1:end) OFDM_frame_delay];


        %%%Receiver 

        %Remove CP
        %OFDM_frame_semcp = OFDM_frame_cp(NCP+1:end);
        OFDM_frame_semcp = OFDM_frame_cp2(NCP+1:end); %remove CP

        %Doppler effect
           for k=1:N
              Frx_new(k,l)=OFDM_frame_semcp(k)*exp(1i*2*pi*To*fd(h)*(l-1));
           end
         
    end
    
    Frx=Frx+Frx_new;

 
%FFT for every column of the matrix
Frx_final=fft(Frx);


%Reshape 1xN*M to [N,M]
Ftx=Ftx(2:end); %delete first element ('0');
Ftx_final=reshape(Ftx,[N,M]); %FTx



%Transfer function F

for n=1:N
   for m=1:M
       F(n,m)=Frx_final(n,m)/Ftx_final(n,m);
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
hold on

end

figure(2);
imagesc(x,y,Per)
title('Periodogram')
xlabel('Relative speed m/s')
ylabel('Distance m')
colorbar


end