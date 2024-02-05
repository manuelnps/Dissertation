
clear all
clc

fc=24e9;  %carrier frequency
lambda=3e8/fc;
B=93.1e6; %signal bandwidth
M=256; %Number of OFDM symbols
N=1024;  %Number of sub-carriers

%M=4; %Number of OFDM symbols
%N=64;  %Number of sub-carriers

% Sub-carrier separation
Tsymbol = N/B ; % Useful symbol duration
Tcp=(1/4)*Tsymbol;
To=Tsymbol+Tcp;

spacing=1/Tsymbol;

Ts = 1/(spacing*N); %Sampling time
fs = 1/Ts; %frequency sampling
NCP = round(Tcp/Ts);

data_mod2=zeros(1,N);

%MIMO PxQ antennas
Q=4; %number of receiver antennas
d_antenna=lambda/2;


%Target Parameters
c0=3e8;
%Resolution and Unambiguous Range/Velocity
resolution_R=c0/(2*spacing*N); %Resolution Range
resolution_v=c0/(2*fc*Tsymbol*M);%Resolution Velocity
resolution_a = lambda/(Q*(0.5*lambda));

%H>1 
%R= m , v= m/s
d=[10]; %Distance m

mind = min(d)-5;
maxd = max(d)+5;

v=[23]; %velocity m/s

minv = min(v)-5;
maxv = max(v)+5;
%fd=[0 0 0 0];
Angle=[70]; %Degrees

mina = -90;
maxa = 90;

angle_r = [-90:90];

for h=1:length(d) %for h targets we'll have 4 different delays because we have for differently spaced receiving antennas (Q)
    for q=1:Q 
         delay(h,q)=2*d(h)/c0+(q-1)*lambda*sind(Angle(h))/(2*c0); %relative delay added and changed it for Angle (3.16 pag 52 and 4.10
         % page 64)
         fd(h,q)=(2*(fc/c0))*v(h); %Doppler Shift
    end
end

%Periodogram
Nper=8*N;
Mper=8*M;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for h=1:length(d) 
    
    for q=1:Q
            
            Ftx=zeros(1,1);
            s_ifft_comp=zeros(1,N);
            F=zeros(N,M);
            Frx_final=zeros(N,M);
            Frx=zeros(N,M);
            Frx_new=zeros(N,M);

                for l=1:M
    
                %%%%Transmitter
    
                    %Data
                    nbits=N*log2(2); %BPSK
                    data=round(rand(1,nbits)); %data generation
                    data_mod = -data*2+1; %data modulation
                    Ftx=[Ftx data_mod];
                    
                    %% -> ciclo for
                    %OFDM frame + CP
                    data_mod2=data_mod.*(exp(-2*j*pi*delay(h,q)*spacing*(0:(N-1)))*exp(-2*1j*pi*fc*delay(h,q)));  %(3.14 Page 52.)
                    
                    
    
                    % Inverse FOURIER Transform
                    s_ifft_comp=ifft(data_mod2);
    
                    % Add CP
    
                    OFDM_frame_cp2 = [s_ifft_comp(end-NCP+1:end) s_ifft_comp];
                 
                %%%%Receiver 
    
                    %Remove CP
                    OFDM_frame_semcp = OFDM_frame_cp2(NCP+1:end); %remove CP
    
                    %Doppler effect
                       for k=1:N
                          Frx(k,l)=OFDM_frame_semcp(k)*exp(1i*2*pi*To*fd(h,q)*(l-1)); %(4.16 page 66)
                       end
    
                end
    
    
            %FFT 
            Frx_final=fft(Frx); %
    
            %Reshape 1xN*M to [N,M]
            Ftx=Ftx(2:end); %delete first element ('0');
            Ftx_final=reshape(Ftx,[N,M]); %FTx
    
             %As we put the data in front of the last data 
             %we need to put it back in matrix form 
    
            %Transfer function F
            for n=1:N
               for m=1:M
                   F_new(n,m)=Frx_final(n,m)/Ftx_final(n,m); %--> Gq(k,l)=fft()*ifft() [NxM]
               end
            end

            %f_new = frx./ftx
            
            F=F+F_new; %target
            
        %end
        
        G(:,:,q)=F_new; %Saving the contents of G(k,l)- k is the index of the subcarrier, and l the ofdm symbol- on a 3D matrix
                        %with q,being the index of the antenna.
         
    end
    
    %%
    % Two-dimension Periodogram
    Per=zeros(Nper,Mper);
    Aux=zeros(N,Mper);
    
    
    Aux=fft(F.',Mper).';
    Per=ifft(Aux,Nper);
    Per=fftshift(Per);
    Per=(1/(N*M)).*abs(Per).^2;
    
    %Plot
    y=[(-Nper/2):((Nper/2)-1)].*(c0/(2*spacing*Nper)); %(3.20, page 53)
    x=[(-Mper/2):((Mper/2)-1)].*(c0/(2*fc*To*Mper)); %(3.25, page 53)
    
    figure(1)
    contour(x,y,Per) 
    title('Periodogram')
    xlabel('Relative speed m/s')
    ylabel('Distance m')
    axis([minv maxv mind maxd])
    colorbar
    hold on;


     for l=1:M
        for k=1:N
            for a=1:180
                %%criar array -90:+90 /somar+91
                G_new=[G(k,l,1);G(k,l,2);G(k,l,3);G(k,l,4)];
                B=exp(2*j*pi*sind((a-1))*(0:(Q-1))*(d_antenna/lambda))'; %the index is the angle "a" not Angle                
                % H(k,l,angle)=sum(G_new(:).*B(:));
                H(k,l,a)=dot(B,G_new);  %the dot product conjugates the first, not the second term


            end
        end
    end
  
        
    distAng = squeeze(H(:,1,:)); %To calculate the 2D periodogram, i used squeeze to remove the dimension
                                 %correspondant to the speed.
    Aper=8*Angle(h);

    Per2 = zeros(Nper, Aper);
    Aux2 = zeros(N,Aper);
    
     %(Eq 5.7, page 69)
    Per2=ifft(distAng,Nper); %As there is no need to do the fft, we can suppress the first lines, and do fftshift with dim = 1.
    Per2=fftshift(Per2,1);
    Per2=(1/(N)).*abs(Per2).^2;
    
    y2=[(-Nper/2):((Nper/2)-1)].*(c0/(2*spacing*Nper)); %(3.20, page 53)
    x2=[0:179]; 
    
    % % Plot Periodogram with Angle and Distance Information
    figure(2);
    contour(x2,y2,Per2) 
    title('Periodogram')
    xlabel('Angle º')
    ylabel('Distance m')
    % Adjust the axis limits as needed
    axis([0 360 mind maxd])
    colorbar;

    hold on

end
    hold off;
    hold off

