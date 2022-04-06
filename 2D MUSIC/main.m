clc; 
clear all;
%--------------------------------------------- 
%--------OFDM Parameters--- 
fc = 2.412 * 10 ^ 9;
%BPSK Modulation
Ntot = 64; %FFT size or total number of subcarriers (used + unused) 64 
N_s_data = 48; %Number of data subcarriers  
N_s_pilots = 4 ; %Number of pilot subcarriers  
ofdmBW = 20 * 10 ^ 6 ; % OFDM bandwidth 
%--------Derived Parameters-------------------- 
deltaF = ofdmBW/ Ntot; %= 20 MHz/ 64 = 0.3125 MHz 
Tfft = 1/ deltaF; % IFFT/ FFT period = 3.2us 
Tgi = Tfft/ 4;% Guard interval duration and also the duration of cyclic prefix
Tsignal = Tgi + Tfft; %duration of BPSK-OFDM symbol 
Ncp = Ntot* Tgi/ Tfft; %Number of symbols allocated to cyclic prefix 
Nst = N_s_data + N_s_pilots; %Number of total used subcarriers 
nBitsPerSym = Nst; %For BPSK the number of Bits per Symbol is same as num of subcarriers 
%-----------------Transmitter-------------------- 
s = 2*randi([0 1], 1, Nst)-1; %Generating Random Data with BPSK modulation 
%IFFT block 
%Assigning subcarriers from 1 to 26 (mapped to 1-26 of IFFT input) 
%and -26 to -1 (mapped to 38 to 63 of IFFT input); 
%Nulls from 27 to 37 and at 0 position 
X_Freq =[ zeros( 1,1) s( 1: Nst/ 2) zeros( 1,11) s( Nst/ 2 + 1: end)]; 
% Assuming that the data is in frequency domain and converting to time domain 
% and scaling the amplitude appropriately
x_Time = Ntot/ sqrt( Nst)* ifft( X_Freq); 
%Adding Cyclic Prefix 
ofdm_signal =[ x_Time( Ntot-Ncp + 1: Ntot) x_Time]; %Generation of the OFDM baseband signal complete
Tsym = Tsignal/(Ntot+Ncp); % duration of each symbol on the OFDM signal
t=Tsym/50:Tsym/50:Tsym; % define a time vector

Carr = BPConv(ofdm_signal,t,fc);

y_tx = Carr;
%NCarr = circshift(y_tx,delay);
%NCarr = delayseq(y_tx',delay)';
rx_carr = y_tx;
c = physconst('LightSpeed');
lam = c/fc;
%  Element Spacing
d = lam;
%  Number of Elements   
N = 4;
theta = [35];
beta = 2*pi;
%beta=2*pi/wavelength; 
phi=beta*(d/lam)*sin(theta*pi/180);
M = length(theta);
%  Steering Vectors
for i=1:M
    for k=1:N
        SteeringVector(k,i)= exp((k-1)*1i*phi(i));
    end
end
rx_carr = SteeringVector.*rx_carr;
%------Channel------
snr = 15;
rx_carr = awgn(rx_carr,snr, 'measured');
delay = 1e-7;
%%-----------------Receiver---------------------- 
%I-Q or vector down-conversion to recover the OFDM baseband signal from the
%modulated RF carrier
lengthofdm = length(ofdm_signal);
r = BBConv(rx_carr,fc,t,N,lengthofdm,Tsym);

for k = 1:N
           for l = 1:26
                arr(l,k) = exp(2*pi*-1i*l*deltaF*delay);
           end
           for l = 38:64
                arr(l-12,k) = exp(2*pi*-1i*l*deltaF*delay);
           end
end
r1 = reshape(r,[80,N]).';
%Removing cyclic prefix 
r_Parallel1 = r1(:,(Ncp + 1:(Ntot + Ncp))); 
%FFT Block 
for i = 1:N
    r_Time(i,:) = sqrt(Nst)/ Ntot*(fft(r_Parallel1(i,:))); 
end

%r_Time = sqrt(Nst)/ 64*(fft(r_Parallel(2,:)));
%Extracting the data carriers from the FFT output 
R_Freq1 = r_Time(:,[( 2: Nst/ 2 + 1) (Nst/ 2 + 13: Nst + 12)]).';
R_Freq1 = arr.*R_Freq1;
%s1 = circshift(s,1);
for i = 1:N
    CSI(:,i)= R_Freq1(:,i)./s.';
end
R_Freq = reshape(R_Freq1,1,[]);
CSI = reshape(CSI,1,[]);

%-------------------------------------------- 
[degrees,delay] = MUSIC(CSI,M,N,d,lam,deltaF);
