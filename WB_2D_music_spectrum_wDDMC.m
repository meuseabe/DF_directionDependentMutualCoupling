% Bu iþlev iki boyutta (yanca , yükseliþ) music spectrumu hesaplar.
% Ahmet M. Elbir. 04.2014
function phi=WB_2D_music_spectrum_wDDMC(Rt,n,k0,r,L,aci_bas,aci_bit,step,sensor_pozisyon,el_bas,el_bit,el_step,AOA)
cj = sqrt(-1);
el_index2 = (el_bas+el_step): el_step: el_bit;
K = size(el_index2,2);
phi = zeros(L,K);
for i = 1%:size(AOA,1)
    Rt = AOA(i,1).R; % covariance matrix for each frequency bin.
%     C = AOA(i,1).C;
    k0 = 2*pi/AOA(i,1).lam; % wavelength for each frequency bin.
    [m,N]=size(Rt);
    % gürültü alt uzayý öz vektörler bulunur.
    [U,D,V]=svd(Rt);
    G=U(:,n+1:m);
    GG = G*G';
    el_index=1;
    
    for k = el_index2
        index = 1;
        for l = (aci_bas+step):step:(aci_bit)
            
            aci=l*pi/180; % yanca.
            el=(k)*pi/180; % yükseliþ.
            
            for kk =1:n
                if l >= AOA(1,1).Sector(1,1) && l <= AOA(1,1).Sector(1,2)
                    C = AOA(1,1).C;
                elseif l >= AOA(2,1).Sector(1,1) && l <= AOA(2,1).Sector(1,2)
                    C = AOA(2,1).C;
                else
                    C = eye(m);
                end
            end
            
            kaynak = k0*[sin(el)*cos(aci.')   sin(el)*sin(aci.') cos(el)];
            a = [ exp(cj*sensor_pozisyon*kaynak.')];
            phi0(index,el_index)=(1/abs(a'*C'*GG*C*a));
            index = index + 1;
        end
        el_index= el_index+1;
    end
    phi0 = double(phi0); % pseudo-spectrum.

    phi = phi + phi0;
end