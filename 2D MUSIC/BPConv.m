function [Carr,t] = BPConv(ofdm_signal,t,fc)
%%----------------------Up-conversion to RF----------------------------
    Carr_real=[];
    Carr_imag=[];
    Carr = [];
    for n=1:length(ofdm_signal)
        Carr_real = real(ofdm_signal(n))*cos(2*pi*fc*t); %modulate the real part on a cosine carrier
        Carr_imag = imag(ofdm_signal(n))*sin(2*pi*fc*t); %modulate the imaginary part on a sine carrier
     %   Carr_real=[Carr_real Carr_real];
     %   Carr_imag=[Carr_imag Carr_imag];
        Carr = [Carr Carr_real+Carr_imag];
    end

end
