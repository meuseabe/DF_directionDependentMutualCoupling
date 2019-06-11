function [AF AN ] = generate_test_data(near_az_index,near_el_index,...
    range,M,lamda,POS,far_az_index,far_el_index,NS)
cj = sqrt(-1);
for ni = 1: length(near_az_index)
    azn = near_az_index(ni);
    eln = near_el_index(ni);
    R = range(ni);
    K_n=[sind(eln)*cosd(azn)   sind(eln)*sind(azn) cosd(eln)];
    for i = 1: M
        an_test(i,ni) = exp(-cj* (2*pi/lamda*R)*(sqrt(1- (2/R)*POS...
            (i,:)* K_n' + sum(POS(i,:).^2)/R^2 )-1));
    end
end

for i = 1: length(far_az_index)
    K_n=[sind(far_el_index(i))*cosd(far_az_index(i))   ....
        sind(far_el_index(i))*sind(far_az_index(i)) cosd(far_el_index(i))];
    af_test(:,i) = exp (cj* (2*pi/lamda*POS*K_n'));
    %     af_test(:,i) = af_test(:,i)/af_test(1,i);
end



AF = af_test(:,1:NS(1)); AN = an_test(:,1:NS(2));