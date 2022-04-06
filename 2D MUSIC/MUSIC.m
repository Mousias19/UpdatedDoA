function [degrees,delay] = MUSIC(CSI,M,N,d,lam,deltaF)

%   MUSIC
%   Auto-correlation matrix
    Rxx = CSI'*CSI/52;
%   eigenvalue decomposition
    [Vi,Li] = eig(Rxx);
%   sort in descending order
    [L,I] = sort(diag(Li),'descend');
    V = Vi(:,I);
%   Signal Subspace
    Ps = V(:,1:M)*V(:,1:M)';
%   Noise Subspace
    Pn = V(:,1+M:N)*V(:,1+M:N)';
    theta1=[0:90];
    tau=[0:1e-9:1e-6];
    for i=1:length(theta1)
        phi1=2*pi*(d/lam)*sin(theta1(i)*pi/180);
        B =zeros([N 1]);
        for k=1:N
            B(k,1)= (exp((k-1)*1i*phi1));
        end
        B = B.';
       for j=1:length(tau)
           for l = 1:26
                A(:,l) = exp(2*pi*-1i*l*deltaF*tau(j));
           end
           for l = 38:64
                A(:,l-12) = exp(2*pi*-1i*l*deltaF*tau(j));
           end
            A1 = kron(B,A)';
            PMUSIC(i,j)= 1./(sum(abs(A1'*Pn*A1)));
       end
    end
    

    [resSorted, orgInd] = sort(PMUSIC,'descend');
    DOAs = orgInd(1:M,1);
    
    [C,I] = max(PMUSIC(:));
    [I1,I2] = ind2sub(size(PMUSIC),I);
    deg = I1-1;
    del = I2-1;
                                      
    figure(1);
    [X,Y] = meshgrid(theta1,tau);
    surf(X,Y,(PMUSIC)')
    shading interp 
    colorbar    
    xlabel('Angle [degrees]');
    ylabel('Delay (s)');
    zlabel('PMUSIC');
    
    degrees = deg;
    delay = del;
end
