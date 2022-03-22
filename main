clc;
clear all;

numFrames = 1; % 10 ms frames
snr = 70; % SNR in dB
fc = 1e8;
LightSpeed = physconst('LightSpeed');   
 %  lambda (λ)
wavelength = LightSpeed/fc; 
 %  Element Spacing
d = 0.5*wavelength;
 %  Number of Elements   
M = 4;
 %  Number of Sources
N = 1;

% Create UE/carrier configuration
ue = nrCarrierConfig;
ue.NSizeGrid = 52;          % Bandwidth in number of resource blocks (52RBs at 15kHz SCS for 10MHz BW)
ue.SubcarrierSpacing = 15;  % 15, 30, 60, 120, 240 (kHz)
ue.CyclicPrefix = 'Normal'; % 'Normal' or 'Extended'

nTxAnts = 1;  % Number of transmit antennas (1,2,4)
nRxAnts = 1;  % Number of receive antennas
nLayers = min(nTxAnts,nRxAnts);

% Configure a periodic multi-port SRS and enable frequency hopping
srs = nrSRSConfig;
srs.NumSRSSymbols = 4;          % Number of OFDM symbols allocated per slot (1,2,4)
srs.SymbolStart = 8;            % Starting OFDM symbol within a slot
srs.NumSRSPorts = nTxAnts;      % Number of SRS antenna ports (1,2,4).
srs.FrequencyStart = 0;         % Frequency position of the SRS in BWP in RBs
srs.NRRC = 0;                   % Additional offset from FreqStart specified in blocks of 4 PRBs (0...67)
srs.CSRS = 14;                  % Bandwidth configuration C_SRS (0...63). It controls the allocated bandwidth to the SRS
srs.BSRS = 0;                   % Bandwidth configuration B_SRS (0...3). It controls the allocated bandwidth to the SRS
srs.BHop = 0;                   % Frequency hopping configuration (0...3). Set BHop < BSRS to enable frequency hopping
srs.KTC = 2;                    % Comb number (2,4). Frequency density in subcarriers
srs.Repetition = 2;             % Repetition (1,2,4). It disables frequency hopping in blocks of |Repetition| symbols
srs.SRSPeriod = [2 0];          % Periodicity and offset in slots. SRSPeriod(2) must be < SRSPeriod(1)
srs.ResourceType = 'periodic';  % Resource type ('periodic', 'semi-persistent','aperiodic'). Use 'aperiodic' to disable inter-slot frequency hopping

csiSubbandSize = 4; % Number of RBs per subband


% Number of slots to simulate
numSlots = numFrames*ue.SlotsPerFrame;

% Total number of subcarriers and symbols per slot
K = ue.NSizeGrid * 12;
L = ue.SymbolsPerSlot;

for nSlot = 0:numSlots-1

    % Update slot counter
    ue.NSlot = nSlot;

    % Generate SRS and map to slot grid
    [srsIndices,srsIndInfo] = nrSRSIndices(ue,srs);
    srsSymbols = nrSRS(ue,srs);

    % Create a slot-wise resource grid empty grid and map SRS symbols
    txGrid = nrResourceGrid(ue,nTxAnts);
    txGrid(srsIndices) = srsSymbols;

    % Determine if the slot contains SRS
    isSRSSlot= ~isempty(srsSymbols);

    % OFDM Modulation
    [txWaveform,waveformInfo] = nrOFDMModulate(ue,txGrid,'CarrierFrequency',fc);
    SampleRate = waveformInfo.SampleRate;
    SamplePeriod = 1/SampleRate;
 %    antenna initialization
 %    antenna = phased.IsotropicAntennaElement('FrequencyRange',[2e8 3e8]);
 %    array = phased.ULA('Element',antenna,'ElementSpacing',d,'NumElements',N);
 %    ang = [17; -20];
 %    doa = [17;60];
    doa = 17;
    %delay = 7e7; %delay
    delay = 3e-6;
    tau = ceil(SampleRate*delay);
    txWaveform = [zeros(tau,N); txWaveform];

%txWaveform 
for k=1:N
    ksi = ceil([0:M-1]'*d*cosd(doa)*(1/LightSpeed)*SampleRate);
end

maxksi = max(ksi);

for h = 1:M
 x(:,h) = [zeros(ksi(h),1); txWaveform; zeros(maxksi-ksi(h),1);];
end
 %x(1:15360,:) = txWaveform1(1:15360,:);

 %  Steering Vectors
 %       for k=1:N
 %           SteeringVector(:, k) = exp(-1i*2*pi*fc*d*(cosd(doa(k))*(1/LightSpeed)*[0:M-1]'+delay)); 
 %       end

 %  Wave arriving at the specified angle
 %   x = SteeringVector*txWaveform.';
 %   x = x.';
 %  OFDM Demodulation
    rxGrid = nrOFDMDemodulate(ue,x,'CarrierFrequency',fc);
 %  Check Normal/Extended Cyclic Prefix
    if strcmp(ue.CyclicPrefix,'normal')
        %   Repeat for SRS Number in a single slot
        for i=9:9+srs.NumSRSSymbols-1
            rxGridEst=zeros(K,M);
            C = permute(rxGrid,[1 3 2]);
            rxGridEst(:,:) = C(:,:,i);
            rxGridEst = rxGridEst.';
            [degrees,delay] = MUSIC(rxGridEst,M,d,N,fc,K,LightSpeed);
        end 
    elseif strcmp(ue.CyclicPrefix,'Extended')
        for i=6:6+srs.NumSRSSymbols-1
            rxGridEst=zeros(K,N);
            C = permute(rxGrid,[1 3 2]);
            rxGridEst(:,:) = C(:,:,i);
            A=1;
            sigma = sqrt((A^2)/(2*10^(snr/10)));
            B = (sigma^2)*(randn(size(rxGridEst))+1i*randn(size(rxGridEst)))/sqrt(2);
            rxGridEst = rxGridEst + B;
            rxGridEst = rxGridEst.';
            [degrees] = MUSIC(rxGridEst,M,d,N,fc,K,LightSpeed);
        end
    else 
        disp('Unexpected Cyclic Prefix Value')
    end
    
end
