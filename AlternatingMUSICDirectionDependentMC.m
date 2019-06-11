function [outputArg] = AlternatingMUSICDirectionDependentMC(inputArg)


%% initials
c = inputArg.c;
Rt = inputArg.Rt;
U = inputArg.U;
Z = inputArg.Z;
i_max = inputArg.i_max;
lamda = inputArg.lamda;
far_az = inputArg.far_az;
far_el = inputArg.far_el;
POS = inputArg.POS;
AOA = inputArg.AOA;
%%
K = length(far_az);
NS = [K 0];
M = size(Rt,1);
Mb = length(c);
k0 = 2*pi/lamda;
near_az_index = far_az;
near_el_index = far_el;
for k = 1:K
    C_est(k,1:M,1:M) = eye(M);
end
for i =1:i_max
    AOA(1,1).R = Rt;
    AOA(1,1).lam = lamda;
    for k = 1:K
        AOA(k,1).C = squeeze(C_est(k,:,:));
        AOA(k,1).Sector = [far_az(k) - 5 far_az(k) + 5];
    end
    %AOA(1,1).el_range = [min(far_el) - 10 max(far_el) + 10];
    [result_MUSIC_w_MC A_w_MC] = WB_2D_MUSIC_MULTIPLE_wDDMC(Rt,k0,K,POS,AOA);
    for k = 1:K
        [val1,ind1] = min(abs(result_MUSIC_w_MC(1,:) - far_az(k)));
        param.it(i).estimated_az(1,k) = result_MUSIC_w_MC(1,ind1);
        param.it(i).estimated_el(1,k) = result_MUSIC_w_MC(2,ind1);
        param.it(i).estimated_A(:,k) = A_w_MC(:,ind1);
    end
    [A_est AN ] = generate_test_data(near_az_index,near_el_index,...
        ones(1,K),M,lamda,POS,param.it(i).estimated_az(1,:),param.it(i).estimated_el(1,:),NS);
    error_az(i,:) = (param.it(i).estimated_az - far_az);
    error_el(i,:) = (param.it(i).estimated_el - far_el);
    absSumError(i,:) = [norm(error_az(i,:))];
    [c_est C_est] = estimate_MCM_DirectionDependent(ones(M,1),A_est,U,Z);
    param.it(i).estimated_c = c_est;
    param.it(i).estimated_C = C_est;
    error_MCM(i,1) = sqrt(norm( c(:) - c_est,2)^2/Mb);
    error = [ error_az error_el error_MCM  ]
    
    errorPlotDOA(i,2) = norm(error_az(i,:))/sqrt(2);
    %errorPlotMC = abs(error_MCM);
    %% error percentage.
    for k = 1:K
        errorPercent0(k,1) = abs(param.it(i).estimated_az(k) - far_az(k) );
        errorPercent0(k,2) = abs(param.it(i).estimated_az(k) - far_az(k) );
        errorPercent1(k,:) = errorPercent0(k,:)/far_az(k)*100;
    end
    errorPercentDOA(i,1) = mean(errorPercent1(:,1));
    errorPercentDOA(i,2) = mean(errorPercent1(:,2));
    errorPlotDOA = errorPercentDOA(:,1);
    
    
    for m = [ 2:Mb Mb+2:Mb+Mb]
        errorPercentMC0Real(m,1) = abs(real(param.it(i).estimated_c(m)) - real(c(m)) );
        errorPercentMC0Real(m,2) = abs(real(param.it(i).estimated_c(m)) - real(c(m)) );
        errorPercentMC1Real(m,:) = errorPercentMC0Real(m,:)/abs(real(c(m)))*100;
        
        errorPercentMC0Imag(m,1) = abs(imag(param.it(i).estimated_c(m)) - imag(c(m)) );
        errorPercentMC0Imag(m,2) = abs(imag(param.it(i).estimated_c(m)) - imag(c(m)) );
        errorPercentMC1Imag(m,:) = errorPercentMC0Imag(m,:)/abs(imag(c(m)))*100;
    end
    
    errorPlotMC(i,1) = mean([errorPercentMC1Real([ 2:Mb Mb+2:Mb+Mb],1);errorPercentMC1Imag([ 2:Mb Mb+2:Mb+Mb],1)]);
   %% 
    
    %[c c_est]
%     if i >2 && norm(error_az(i,:) - error_az(i-1,:) ,2)   <0.001 ...
%             && norm(error_az(i-1,:) - error_az(i-2,:) ,2)   <0.001
%         break;end
end
%% iteration figure
    figure(10)
    subplot(211)
    plot(1:i,errorPlotDOA,'LineWidth',2,'Marker','v')
    xlabel('ITERATION NUMBER')
    ylabel('ERROR PERCENT, DOA')
    title('SNR=12dB, T=200')
    grid on
    axis tight
    subplot(212)
    plot(1:i,errorPlotMC,'LineWidth',2,'Marker','v')
    xlabel('ITERATION NUMBER')
    ylabel('ERROR PERCENT, MC')
    axis tight
    grid on
%% outputs
if inputArg.da == 1
[minVal, minInd] = min(absSumError);
else
minInd =i;
end
outputArg.estimated_az = param.it(minInd).estimated_az;
outputArg.estimated_el = param.it(minInd).estimated_el;
outputArg.estimated_C = param.it(minInd).estimated_C;
outputArg.estimated_c =  param.it(minInd).estimated_c;
