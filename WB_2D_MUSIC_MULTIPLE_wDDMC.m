% bu iþlev verilen kovaryan matrisi için iki boyutta (yanca,yükseliþ)
% çoklu kaynak durumu için açý kestirimi yapar.
% Ahmet M. Elbir. 04.2014
function [result_MUSIC_w_MC A phi_b xx]=WB_2D_MUSIC_MULTIPLE_wDDMC(Rt,k0,kaynak_sayisi,UCA_pozisyon,AOA)

n=kaynak_sayisi;
% burada etkili bir yanca ve yükseliþ taramasý için sýnýrlar ve tarama aralýðý belirlenir.
aci_bas=-2;aci_bit=180+2;
L = 64*(aci_bit-aci_bas);step=(aci_bit-aci_bas)/L; % yanca için.

el_bas=AOA(1,1).el_range(1,1);el_bit=AOA(1,1).el_range(1,2);
K = round(el_bit-el_bas);el_step=(el_bit-el_bas)/K; % yükseliþ için.
% kabaca spectrum taramasý
phi_b = WB_2D_music_spectrum_wDDMC(Rt,n,k0,1,L,aci_bas,aci_bit,step,UCA_pozisyon,el_bas,el_bit,el_step,AOA);
% spectrumu çizdir.
% figure(112)
xx = linspace(aci_bas,aci_bit,L);
% yy = linspace(el_bas,el_bit,K);
% [XX,YY] = meshgrid(xx,yy);
% mesh(XX,YY,phi_b'/max(max(phi_b)))
% xlabel('AZIMUTH, (Deg.)')
% ylabel('ELEVATION, (Deg.)')
% zlabel('MUSIC PSEUDO-SPECTRUM')
% set(gca,'zscale','log');
% rotate3d on
% hold on
% axis tight
% view(76,27)
% view(0,90)
% hold off
% spectrum için tepe noktalarýný bul.
[peak_m,peak_n]=findpeaks2d(phi_b,kaynak_sayisi);

for kk = 1:size(peak_m,2)
peak_most(kk) = phi_b(peak_m(kk),peak_n(kk));
az_index_peaks(kk) = peak_m(kk);
el_index_peaks(kk) = peak_n(kk);
end
% peak_most
[peak_found index_found] = sort(peak_most,'descend');
if kaynak_sayisi>size(peak_m,2)
    n = size(peak_m,2); % bu durumda kaynaklar fully coherent. kaynak sayisi kadar tepe yok.
else
    n = size(peak_m,2); % diðer durumda da tüm tepe noktalarýný raporla.
end
for r = 1:n
    az_index (r) = az_index_peaks(index_found(r));
    el_index (r) = el_index_peaks(index_found(r));
UCA_MUSIC(r,:) = [(aci_bas + az_index(r)*step) ;(el_bas + el_index(r)*el_step)];
end
% detaylý arama yap
for r = 1:n
    step=.5;
    el_step=.5;
    L = 50;aci_bas(r)=UCA_MUSIC(r,1)-step;aci_bit(r)=UCA_MUSIC(r,1)+step;step=(aci_bit(r)-aci_bas(r))/L;
    K = 50; el_bas(r)=UCA_MUSIC(r,2)-el_step;el_bit(r)=UCA_MUSIC(r,2)+el_step;el_step=(el_bit(r)-el_bas(r))/K;
    % detaylý spectrum taramasý yap.
    phi_2 = WB_2D_music_spectrum_wDDMC(Rt,kaynak_sayisi,k0,1,L,aci_bas(r),aci_bit(r),step,UCA_pozisyon,el_bas(r),el_bit(r),el_step,AOA);
    [max_value_x,max_index_x]=max(phi_2);
    [max_value_y,max_index_y]=max(max(phi_2));
    % tepe deðerleri al.
    UCA_MUSIC(r,:)=[aci_bas(r) + max_index_x(1)*step, el_bas(r)+el_step*max_index_y];
end

for r = 1:n
    step=.015;
    el_step=0.015;
    L = 50;aci_bas(r)=UCA_MUSIC(r,1)-step;aci_bit(r)=UCA_MUSIC(r,1)+step;step=(aci_bit(r)-aci_bas(r))/L;
    K = 50; el_bas(r)=UCA_MUSIC(r,2)-el_step;el_bit(r)=UCA_MUSIC(r,2)+el_step;el_step=(el_bit(r)-el_bas(r))/K;
    % detaylý spectrum taramasý yap.
    phi_2 = WB_2D_music_spectrum_wDDMC(Rt,kaynak_sayisi,k0,1,L,aci_bas(r),aci_bit(r),step,UCA_pozisyon,el_bas(r),el_bit(r),el_step,AOA);
    [max_value_x,max_index_x]=max(phi_2);
    [max_value_y,max_index_y]=max(max(phi_2));
    % tepe deðerleri al.
    UCA_MUSIC(r,:)=[aci_bas(r) + max_index_x(1)*step, el_bas(r)+el_step*max_index_y];
end

for r = 1:n
    step=0.0015;
    el_step=0.0015;
    L = 50;aci_bas(r)=UCA_MUSIC(r,1)-step;aci_bit(r)=UCA_MUSIC(r,1)+step;step=(aci_bit(r)-aci_bas(r))/L;
    K = 50; el_bas(r)=UCA_MUSIC(r,2)-el_step;el_bit(r)=UCA_MUSIC(r,2)+el_step;el_step=(el_bit(r)-el_bas(r))/K;
    % detaylý spectrum taramasý yap.
    phi_2 = WB_2D_music_spectrum_wDDMC(Rt,kaynak_sayisi,k0,1,L,aci_bas(r),aci_bit(r),step,UCA_pozisyon,el_bas(r),el_bit(r),el_step,AOA);
    [max_value_x,max_index_x]=max(phi_2);
    [max_value_y,max_index_y]=max(max(phi_2));
    % tepe deðerleri al.
    UCA_MUSIC(r,:)=[aci_bas(r) + max_index_x(1)*step, el_bas(r)+el_step*max_index_y];
end

% for r = 1:n
%     step=0.0015;
%     el_step=0.0015;
%     L = 50;aci_bas(r)=UCA_MUSIC(r,1)-step;aci_bit(r)=UCA_MUSIC(r,1)+step;step=(aci_bit(r)-aci_bas(r))/L;
%     K = 50; el_bas(r)=UCA_MUSIC(r,2)-el_step;el_bit(r)=UCA_MUSIC(r,2)+el_step;el_step=(el_bit(r)-el_bas(r))/K;
%     % detaylý spectrum taramasý yap.
%     phi_2 = WB_2D_music_spectrum_wMC(Rt,kaynak_sayisi,k0,1,L,aci_bas(r),aci_bit(r),step,UCA_pozisyon,el_bas(r),el_bit(r),el_step,AOA);
%     [max_value_x,max_index_x]=max(phi_2);
%     [max_value_y,max_index_y]=max(max(phi_2));
%     % tepe deðerleri al.
%     UCA_MUSIC(r,:)=[aci_bas(r) + max_index_x(1)*step, el_bas(r)+el_step*max_index_y];
% end
% UCA_MUSIC(:,2) = UCA_MUSIC(:,2);

cj = sqrt(-1);
K= n;...kaynak_sayisi;
result_MUSIC_w_MC = UCA_MUSIC';
% K1 = size(result_MUSIC_w_MC,2);
% if  K1 < K
%     result_MUSIC_w_MC(:,K1+1:K) = kron(ones(1,K-K1),result_MUSIC_w_MC(:,K1));
% end
% for k = 1:K
%     if result_MUSIC_w_MC(1,k) >=180
%         result_MUSIC_w_MC(1,k) = abs(result_MUSIC_w_MC(1,k)-180);
%     end
% end
for r =1:K
    aci_az(1,r) = result_MUSIC_w_MC(1,r);
    aci_el(1,r) = result_MUSIC_w_MC(2,r);
    kaynak = k0*[sind(aci_el(1,r))*cosd(aci_az(1,r).')   sind(aci_el(1,r))*sind(aci_az(1,r).') cosd(aci_el(1,r))];
    A(:,r) = [ exp(cj*UCA_pozisyon*kaynak.')];
end


