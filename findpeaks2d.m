% Bu kod. 2 boyutlu bir matristeki tepe noktalarını bulur.
% Ahmet M. Elbir. 04.2014
function [peak_m,peak_n,peak_value_m,peak_value_n] = findpeaks2d(P,number_of_peaks)

[M,N] = size(P);

peak_m = [];
peak_n = [];

% for m = 2:M-1
%     for n = 2:N-1
%         
%         if      P(m,n)> P(m-1,n-1) &&...
%                 P(m,n)> P(m-1,n) &&...
%                 P(m,n)> P(m-1,n+1) &&...
%                 P(m,n)> P(m,n-1) &&...
%                 P(m,n)> P(m,n+1) &&...
%                 P(m,n)> P(m+1,n-1) &&...
%                 P(m,n)> P(m+1,n) &&...
%                 P(m,n)> P(m+1,n+1)
%             % peak found
%             peak_m = [ peak_m m ];
%             peak_n = [ peak_n n ];
%         end
%         
%     end
% end


for m = 1:M % az
    Pm(m,1) = sum(P(m,:));
end

% find peaks in az.
peak_m0 = [];
peak_value_m = [];
for m = 2:length(Pm)-1
    if Pm(m,1) > Pm(m-1,1) &&...
            Pm(m,1)> Pm(m+1,1)
        peak_m0 = [peak_m0 m];
        peak_value_m = [peak_value_m Pm(m,1)];
    else
    end
end
% peak_m = sort(fpeaks4(Pm',number_of_peaks+100,length(Pm),1));

for k = 1:length(peak_m0)
Pn = P(peak_m0(1,k),:);
[val,ind0] = sort(Pn,'descend');
[peak_n0(k)] = ind0(1);
[peak_value_n0(k)] = Pn(ind0(1));
end

[val_sorted, ind_sorted] = sort(peak_value_n0,'descend');

peak_value_n = peak_value_n0(ind_sorted);
peak_n = peak_n0(ind_sorted);
peak_m = peak_m0(ind_sorted);
% pm = peak_m(1:number_of_peaks);

% pn = peak_n(1:number_of_peaks);
% peak_m
if isempty(peak_m)
    [val,ind0] = sort(P(:),'descend');
    for k = 1:number_of_peaks
         [peak_m(k),peak_n(k)] = find(P==val(k));
    end
end
pt = 5; % report at most pt result.
if length(peak_m) >=pt
peak_m = peak_m(1:pt);
peak_n = peak_n(1:pt);
end