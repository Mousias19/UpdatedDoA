function r = BBConv(rx_carr,fc,t,N,lengthofdm,Tsym)
    r = [];
    r_real= [];
    r_imag = [];
    for k=1:N
        for n=1:1:lengthofdm
        %%XXXXXX inphase coherent dector XXXXXXX
        Z_in=rx_carr(k,(n-1)*length(t)+1:n*length(t)).*cos(2*pi*fc*(t)); %extract a period of the 
        %signal received and multiply the received signal with the cosine component of the carrier signal

        Z_in_intg=(trapz(t,Z_in(1,:)))*(2/Tsym);% integration using Trapizodial rule 
                                        %over a period of half bit duration
        r_real=Z_in_intg;

        %%XXXXXX Quadrature coherent dector XXXXXX
        Z_qd=rx_carr(k,(n-1)*length(t)+1:n*length(t)).*sin(2*pi*fc*(t));
        %above line indicat multiplication ofreceived & Quadphase carred signal

        Z_qd_intg=(trapz(t,Z_qd(1,:)))*(2/Tsym);%integration using trapizodial rule

        r_imag =  Z_qd_intg;   

        r=[r  r_real+1i*(r_imag)]; % Received Data vector
        end
    end

end
