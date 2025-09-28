%% Correction for all slices and tpn with central slice (2D) phase and/or MUSE
acc_MB = 1;
[Ny_SB,Nx_SB,Nz_SB,Ncoil_SB,Nt_SB] = size(kspace_SB_tpn);

Image_para = zeros(Ny_SB,Nx_SB);
Image_para(Ny_SB/2,Nx_SB/2) = 1;
kspace_para = ifftshift(ifft2(ifftshift(Image_para)));
para_PosNeg = zeros(ns4,ns4);
for k1 = 1 : ns4
    kspace_k1 = zeros(size(kspace_para));
    kspace_k1(k1:ns4:Ny_SB,:) = kspace_para(k1:ns4:Ny_SB,:);
    Image_k1 = fftshift(fft2(fftshift(kspace_k1)));
    para_PosNeg(:,k1) = angle(Image_k1(sum(abs(Image_k1),2)>1e-5,Nx_SB/2));
end
para_index = [nshot:-1:1 ns4:-1:nshot+1];

Pos_k1 = transpose(exp(1i.*para_PosNeg(para_index,1)));
Pos_k2 = transpose(exp(1i.*para_PosNeg(para_index,2)));
Neg_k3 = transpose(exp(1i.*para_PosNeg(para_index,3)));
Neg_k4 = transpose(exp(1i.*para_PosNeg(para_index,4)));

FOV = Ny_SB;

Image_pos_k1 = zeros(Ny_SB,Nx_SB,Ncoil_SB);
Image_pos_k2 = zeros(Ny_SB,Nx_SB,Ncoil_SB);
Image_neg_k3 = zeros(Ny_SB,Nx_SB,Ncoil_SB);
Image_neg_k4 = zeros(Ny_SB,Nx_SB,Ncoil_SB);

CorrImg_SB_2shot_2D = zeros(Ny_SB,Nx_SB,Nz_SB,Nt_SB);

for z = 1 : Nz_SB
    for tpn = 1 : Nt_SB
    if z == 1 && tpn == 1
        disp(['Started: correting SB images (2D, MUSE) at slice #' num2str(z) '/' num2str(Nz_SB) ', at tpn: ' num2str(tpn) '/' num2str(Nt_SB)]);
    end
        for c = 1 : Ncoil_SB
            kspace1 = zeros(Ny_SB,Nx_SB);
            kspace2 = zeros(Ny_SB,Nx_SB);
            kspace3 = zeros(Ny_SB,Nx_SB);
            kspace4 = zeros(Ny_SB,Nx_SB);
            
            kspace1(1:ns4:end,:) = kspace_SB_tpn(1:ns4:end,:,z,c,tpn);
            kspace2(2:ns4:end,:) = kspace_SB_tpn(2:ns4:end,:,z,c,tpn);
            kspace3(3:ns4:end,:) = kspace_SB_tpn(3:ns4:end,:,z,c,tpn);
            kspace4(4:ns4:end,:) = kspace_SB_tpn(4:ns4:end,:,z,c,tpn);
            
            Image_pos_k1(:,:,c) = fftshift(fft2(fftshift(kspace1)));
            Image_pos_k2(:,:,c) = fftshift(fft2(fftshift(kspace2)));
            Image_neg_k3(:,:,c) = fftshift(fft2(fftshift(kspace3)));
            Image_neg_k4(:,:,c) = fftshift(fft2(fftshift(kspace4)));
        end
        
        MUSE_phase = zeros(Ny_SB,Nx_SB); % We did not include MUSE reconstruction here

        % Final correction
        for k1 = 1 : Nx_SB % x direction
            LineCorr_s0 = zeros(Ny_SB,1);

            CSP_matrix =  zeros(ns4*Ncoil_SB*Ny_SB/ns4,Ny_SB*acc_MB); % csp and phase matrix
            Imat = zeros(ns4*Ncoil_SB*Ny_SB/ns4,1); % Odd-even matrix
            
            for k2 = 1 : Ny_SB/ns4  % y direction
                
                p_s0 = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4];
                p_s0(p_s0>FOV) = p_s0(p_s0>FOV) - FOV;                
                
                pc_r = exp(1i.*transpose(PhaseMap_2D_SB(p_s0,k1)));
                
                MUSE_2shot = exp(1i.*transpose(MUSE_phase(p_s0,k1)));

                for c = 1 : Ncoil_SB
                    CSP = transpose(   CSP_SB(p_s0,k1,z,c));
                    CSP_matrix((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4+1,(k2-1)*ns4 + 1:(k2-1)*ns4 + 4) = 1/ns4.*Pos_k1.*pc_r.*CSP;
                    CSP_matrix((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4+2,(k2-1)*ns4 + 1:(k2-1)*ns4 + 4) = 1/ns4.*Pos_k2.*pc_r.*CSP.*MUSE_2shot;
                    
                    CSP_matrix((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4+3,(k2-1)*ns4 + 1:(k2-1)*ns4 + 4) = 1/ns4.*Neg_k3.*CSP;
                    CSP_matrix((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4+4,(k2-1)*ns4 + 1:(k2-1)*ns4 + 4) = 1/ns4.*Neg_k4.*CSP.*MUSE_2shot;
                    
                end
                Imat((k2-1)*Ncoil_SB*ns4 + 1 : ns4 : k2*Ncoil_SB*ns4) = squeeze(Image_pos_k1(k2,k1,:));
                Imat((k2-1)*Ncoil_SB*ns4 + 2 : ns4 : k2*Ncoil_SB*ns4) = squeeze(Image_pos_k2(k2,k1,:));
                Imat((k2-1)*Ncoil_SB*ns4 + 3 : ns4 : k2*Ncoil_SB*ns4) = squeeze(Image_neg_k3(k2,k1,:));
                Imat((k2-1)*Ncoil_SB*ns4 + 4 : ns4 : k2*Ncoil_SB*ns4) = squeeze(Image_neg_k4(k2,k1,:));
            end
            
            Phase = CSP_matrix;
            S = Phase\Imat;
            
            for k2 = 1 : Ny_SB/ns4
                
                p_s0 = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4];
                p_s0(p_s0>FOV) = p_s0(p_s0>FOV) - FOV;
                
                LineCorr_s0(p_s0,1) = S((k2-1)*ns4 + 1:(k2-1)*ns4 + 4,1);
            end
            
            
            CorrImg_SB_2shot_2D(:,k1,z,tpn) = LineCorr_s0;
        end
    end
    if z == Nz_SB && tpn == Nt_SB
        disp(['Completed: correting SB images (2D, MUSE) at slice #' num2str(z) '/' num2str(Nz_SB) ', at tpn: ' num2str(tpn) '/' num2str(Nt_SB)]);
    end
end
