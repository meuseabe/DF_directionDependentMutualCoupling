
clear all

cj = sqrt(-1);
fc = 100e6;
c0 = 300e6;
lamda = c0/fc;
k0 = 2*pi/lamda;
d = .5*lamda;
%%%%%%%%%%%

%% simulation parameters.
M = 9;
num_trial = 10;
SNR_index =20:5:30;
snap_index = [200];
% as = 20;far_az_index = [50 80 80 ] + as;
as = 90;far_az_index = [35 124 60 ];
far_el_index = [90 90 90]-0;
single_DOA = 1;
single_GPM = 0;
single_MCM = 0;
cal_step = .01;
NS = [2 0]; % exact number of FF and NF sources
% Array Type:
% 1: ULA
% 2: UCA
% 3: 2-D RA
% 4: 3-D RA
% 5: Co-prime Array
% 6: Nested Array
Array_Type = 1;
DirectionDependent = 1;
K = sum(NS);
%% generate array positions
UCA_array_spacing = 360 / M;
r = d/(2*sind(UCA_array_spacing/2));
for i = 1:M
    pos_UCA(i,:) = r* [cosd(UCA_array_spacing*(i-1)) sind(UCA_array_spacing*(i-1)) 0];
    pos_ULA(i,:) = [d*(i-1) 10 10];
end
%% Generate Data
for ns = 1:length(snap_index)
    for sc = 1: size(SNR_index,2)
        SNR = SNR_index(sc); % dB
        for nt = 1: num_trial
            T = snap_index(ns); % number of snapshot.

                %%
                POS = pos_ULA;
                M2_upper_limit = M-K-4;
                if single_MCM == 1
                    [c_ULA] = generate_MCM(M2_upper_limit);
                    single_MCM = 2;
                elseif single_MCM == 0
                    [c_ULA] = generate_MCM(M2_upper_limit);
                end
                Mb = M2_upper_limit;
                if DirectionDependent == 1
                    for k =1:K
                        cMc(:,k) = [1; c_ULA(2:end,1)*k*.8];
                        %cMc(:,k) = [1; .2 + 1i*.2; .2 + 1i*.1; .1 + 1i*.1; + .04 + 1i*.03];
                        c_ULA = cMc(:,k);
                        Zc_ULA(1).E = eye(M); C_ULA = eye(M);%c_UCA(1);
                        for m = 2:Mb
                            e1 = zeros(1,M);
                            e1(1,m) = 1; e1(1,m) = 1;
                            Zc_ULA(m).E = toeplitz(e1,e1);
                            clear e1
                            C_ULA = C_ULA + Zc_ULA(m).E*c_ULA(m);
                        end
                        Z = Zc_ULA;
                        Mb = length(Zc_ULA);
                        c(:,k) = c_ULA;
                        C(k,:,:) = C_ULA;
                    end
                else
                end

            %% generate test data
            far_az = far_az_index(1:NS(1)) + 0.00*(rand(1,1)-0.5);
            far_el = far_el_index(1:NS(1)) + 0.00*(rand(1,1)-0.5);

            [AF AN ] = generate_test_data(far_az,far_el,...
                far_az,M,lamda,POS,far_az,far_el,NS);

            AA = sum([AF AN],2);
            %% measurement generation.
            % independent sources.
            S = 0.5 + abs(randn(K,T));

            AS = [AF AN]*S;
            Mismatch_matrix = C;
            Mismatch_vector = c;
            Zmismatch = Z;
            X = zeros(M,T);
            for k = 1:K
                X = X + squeeze(Mismatch_matrix(k,:,:))*AF(:,k)*S(k,:);
            end
            %% save parameters.
            Observation(ns,sc,nt).AS = AS;
            Observation(ns,sc,nt).AF = AF;
            Observation(ns,sc,nt).X = X;
            Observation(ns,sc,nt).B(:,1:T) = awgn(X,SNR);
            Observation(ns,sc,nt).noise(:,1:T) = Observation(ns,sc,nt).B(:,1:T) - X;
            Observation(ns,sc,nt).Zmismatch = Zmismatch;
            Observation(ns,sc,nt).Mismatch_matrix = Mismatch_matrix;
            Observation(ns,sc,nt).Mismatch_vector = Mismatch_vector;
            Observation(ns,sc,nt).far_az = far_az;
            Observation(ns,sc,nt).far_el = far_el;
            Observation(ns,sc,nt).POS = POS;
        end
    end
end



%% Simulations
tic
for ns = 1:length(snap_index)
    T = snap_index(1,ns);
    for sc = 1: size(SNR_index,2)
        count1 = 0; count2 = 0; count3 = 0; count4 = 0; count5 = 0; count6 = 0;
        for nt = 1: num_trial
            %% load data
            %A = Observation(ns,sc,nt).A;
            AS = Observation(ns,sc,nt).AS;
            AF = Observation(ns,sc,nt).AF;
            X = Observation(ns,sc,nt).X;
            B = Observation(ns,sc,nt).B;
            %G = Observation(ns,sc,nt).G;
            %g = Observation(ns,sc,nt).g;
            %D = Observation(ns,sc,nt).D;
            Z = Observation(ns,sc,nt).Zmismatch; Mb = size(Z,2);
            C = Observation(ns,sc,nt).Mismatch_matrix;
            c = Observation(ns,sc,nt).Mismatch_vector;
            noise = Observation(ns,sc,nt).noise;
            far_az = Observation(ns,sc,nt).far_az;
            far_el = Observation(ns,sc,nt).far_el;
            %Af_ID = Observation(ns,sc,nt).I;
            POS = Observation(ns,sc,nt).POS;
            Rt = B*B';
            [U1,V1] = eig(Rt);
            U2 = fliplr(U1);
            U = U2(:,K+1:M);
            
            
            %% Proposed Method, MUSIC + direction dependent MC
            AOA(1,1).R = Rt; % covariance matrix
            AOA(1,1).lam = lamda; % wavelength
            AOA(1,1).el_range = [min(far_el) - 0 max(far_el) + 1];
            i_max = 30;
            inputArgDD.c = c; % MC vector
            inputArgDD.Rt = Rt; % covariance matrix
            inputArgDD.U = U; % noise subspace
            inputArgDD.Z = Z; % to construct C
            inputArgDD.i_max = i_max; % max iter
            inputArgDD.lamda = lamda; % wavelength
            inputArgDD.far_az = far_az; % true az
            inputArgDD.far_el = far_el; % true el
            inputArgDD.POS = POS; % positions of antennas
            inputArgDD.AOA = AOA;
            [outputArgDD] = AlternatingMUSICDirectionDependentMC(inputArgDD);
            RESULTS.AZ1(ns,nt,sc,1:K) = outputArgDD.estimated_az;
            RESULTS.EL1(ns,nt,sc,1:K) = outputArgDD.estimated_el;
            RESULTS.C1(ns,nt,sc,1:K,:,:) = outputArgDD.estimated_C;
            RESULTS.c1(ns,nt,sc,1:Mb*K) = outputArgDD.estimated_c;
            
            
        end
    end
end
