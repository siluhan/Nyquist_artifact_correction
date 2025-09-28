%% Correction for all slices and tpn with central slice phase (2D) and/or MUSE

nshot = 2;
ns4 = nshot*2;
acc_MB = 2;

Image_para = zeros(Ny_MB,Nx_MB);
Image_para(Ny_MB/2,Nx_MB/2) = 1;
kspace_para = ifftshift(ifft2(ifftshift(Image_para)));
para_PosNeg = zeros(ns4,ns4);
for k1 = 1 : ns4
    kspace_k1 = zeros(size(kspace_para));
    kspace_k1(k1:ns4:Ny_MB,:) = kspace_para(k1:ns4:Ny_MB,:);
    Image_k1 = fftshift(fft2(fftshift(kspace_k1)));
    para_PosNeg(:,k1) = angle(Image_k1(sum(abs(Image_k1),2)>1e-5,Nx_MB/2));
end
para_index = [nshot:-1:1 ns4:-1:nshot+1];

Pos_k1 = repmat(transpose(exp(1i.*para_PosNeg(para_index,1))),1,acc_MB);
Pos_k2 = repmat(transpose(exp(1i.*para_PosNeg(para_index,2))),1,acc_MB);
Neg_k3 = repmat(transpose(exp(1i.*para_PosNeg(para_index,3))),1,acc_MB);
Neg_k4 = repmat(transpose(exp(1i.*para_PosNeg(para_index,4))),1,acc_MB);

FOV = Ny_MB;

Image_pos_k1 = zeros(Ny_MB,Nx_MB,Ncoil_MB);
Image_pos_k2 = zeros(Ny_MB,Nx_MB,Ncoil_MB);
Image_neg_k3 = zeros(Ny_MB,Nx_MB,Ncoil_MB);
Image_neg_k4 = zeros(Ny_MB,Nx_MB,Ncoil_MB);

CorrImg_MB_2shot_2D = zeros(Ny_MB,Nx_MB,Nz_MB*acc_MB);
for z = 1 : Nz_MB
    for tpn = 1 : Nt_MB
        if z == 1 && tpn == 1
            disp(['Started: correting MB images (2D phase, MUSE) at slice #' num2str(z) '/' num2str(Nz_MB) ', at tpn: ' num2str(tpn) '/' num2str(Nt_MB)]);
        end
        
        for c = 1 : Ncoil_MB
            kspace1 = zeros(Ny_MB,Nx_MB);
            kspace2 = zeros(Ny_MB,Nx_MB);
            kspace3 = zeros(Ny_MB,Nx_MB);
            kspace4 = zeros(Ny_MB,Nx_MB);
            
            kspace1(1:ns4:end,:) = kspace_MB_tpn(1:ns4:end,:,z,c,tpn);
            kspace2(2:ns4:end,:) = kspace_MB_tpn(2:ns4:end,:,z,c,tpn);
            kspace3(3:ns4:end,:) = kspace_MB_tpn(3:ns4:end,:,z,c,tpn);
            kspace4(4:ns4:end,:) = kspace_MB_tpn(4:ns4:end,:,z,c,tpn);
            
            Image_pos_k1(:,:,c) = fftshift(fft2(fftshift(kspace1)));
            Image_pos_k2(:,:,c) = fftshift(fft2(fftshift(kspace2)));
            Image_neg_k3(:,:,c) = fftshift(fft2(fftshift(kspace3)));
            Image_neg_k4(:,:,c) = fftshift(fft2(fftshift(kspace4)));
        end       

        % We did not include MUSE recon here
        MUSE_phase = zeros(Ny_SB,Nx_SB,2);
        MUSE_phase(:,:,1) = zeros(Ny_SB,Nx_SB);
        MUSE_phase(:,:,2) = zeros(Ny_SB,Nx_SB);        
        
        for k1 = 1 : Nx_MB % x direction
            LineCorr_s0 = zeros(Ny_MB,1);
            LineCorr_s1 = zeros(Ny_MB,1);
            
            CSP_matrix =  zeros(ns4*Ncoil_MB*Ny_MB/ns4,Ny_MB*acc_MB); % csp and phase matrix
            Imat = zeros(ns4*Ncoil_MB*Ny_MB/ns4,1); % Odd-even matrix
            
            for k2 = 1 : Ny_MB/ns4  % y direction
                
                if FOVshift == 1/4
                    p_s1 = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4] - FOV/4;
                    p_s1(p_s1<1) = p_s1(p_s1<1) + FOV;
                    
                    p_s0 = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4];
                    p_s0(p_s0>FOV) = p_s0(p_s0>FOV) - FOV;
                    
                elseif FOVshift == 1/2
                    p_s0 = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4];
                    p_s0(p_s0>FOV) = p_s0(p_s0>FOV) - FOV;
                    
                    p_s1 = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4]+FOV/2;
                    p_s1(p_s1>FOV) = p_s1(p_s1>FOV) - FOV;
                elseif FOVshift == 0
                    p_s1 = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4];
                    p_s1(p_s1>FOV) = p_s1(p_s1>FOV) - FOV;
                    
                    p_s0 = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4];
                    p_s0(p_s0>FOV) = p_s0(p_s0>FOV) - FOV;
                end
                
                pc_r = exp(1i.*transpose([PhaseMap(p_s0,k1,1);PhaseMap(p_s1,k1,2)]));
                
                MUSE_2shot = exp(1i.*transpose([MUSE_phase(p_s0,k1,1);MUSE_phase(p_s1,k1,2)]));
                for c = 1 : Ncoil_MB
                    CSP = transpose([CSP_SB(p_s0,k1,MB_SliceNum(1,z),c);CSP_SB(p_s1,k1,MB_SliceNum(2,z),c)]);
                    CSP_matrix((k2-1)*Ncoil_MB*ns4 + (c-1)*ns4+1,(k2-1)*ns4*2 + 1:(k2-1)*ns4*2 + ns4*2) = 1/ns4.*Pos_k1.*pc_r.*CSP;
                    CSP_matrix((k2-1)*Ncoil_MB*ns4 + (c-1)*ns4+2,(k2-1)*ns4*2 + 1:(k2-1)*ns4*2 + ns4*2) = 1/ns4.*Pos_k2.*pc_r.*CSP.*MUSE_2shot;
                    
                    CSP_matrix((k2-1)*Ncoil_MB*ns4 + (c-1)*ns4+3,(k2-1)*ns4*2 + 1:(k2-1)*ns4*2 + ns4*2) = 1/ns4.*Neg_k3.*CSP;
                    CSP_matrix((k2-1)*Ncoil_MB*ns4 + (c-1)*ns4+4,(k2-1)*ns4*2 + 1:(k2-1)*ns4*2 + ns4*2) = 1/ns4.*Neg_k4.*CSP.*MUSE_2shot;
                    
                end
                Imat((k2-1)*Ncoil_MB*ns4 + 1 : ns4 : k2*Ncoil_MB*ns4) = squeeze(Image_pos_k1(k2,k1,:));
                Imat((k2-1)*Ncoil_MB*ns4 + 2 : ns4 : k2*Ncoil_MB*ns4) = squeeze(Image_pos_k2(k2,k1,:));
                Imat((k2-1)*Ncoil_MB*ns4 + 3 : ns4 : k2*Ncoil_MB*ns4) = squeeze(Image_neg_k3(k2,k1,:));
                Imat((k2-1)*Ncoil_MB*ns4 + 4 : ns4 : k2*Ncoil_MB*ns4) = squeeze(Image_neg_k4(k2,k1,:));
            end
            
            Phase = CSP_matrix;
            S = Phase\Imat;
            
            for k2 = 1 : Ny_MB/ns4
                
                if FOVshift == 1/4
                    p_s1 = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4] - FOV/4;
                    p_s1(p_s1<1) = p_s1(p_s1<1) + FOV;
                    
                    p_s0 = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4];
                    p_s0(p_s0>FOV) = p_s0(p_s0>FOV) - FOV;
                    
                elseif FOVshift == 1/2
                    p_s0 = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4];
                    p_s0(p_s0>FOV) = p_s0(p_s0>FOV) - FOV;
                    
                    p_s1 = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4]+FOV/2;
                    p_s1(p_s1>FOV) = p_s1(p_s1>FOV) - FOV;
                elseif FOVshift == 0
                    p_s1 = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4];
                    p_s1(p_s1>FOV) = p_s1(p_s1>FOV) - FOV;
                    
                    p_s0 = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4];
                    p_s0(p_s0>FOV) = p_s0(p_s0>FOV) - FOV;
                end
                
                LineCorr_s0(p_s0,1) = S((k2-1)*ns4*2 + 1:(k2-1)*ns4*2 + ns4,1);
                LineCorr_s1(p_s1,1) = S((k2-1)*ns4*2 + 5:(k2-1)*ns4*2 + ns4*2,1);
            end
            
            
            CorrImg_MB_2shot_2D(:,k1,MB_SliceNum(1,z),tpn) = LineCorr_s0;
            CorrImg_MB_2shot_2D(:,k1,MB_SliceNum(2,z),tpn) = LineCorr_s1;
        end
    end
    if z == Nz_MB && tpn == Nt_MB
        disp(['Completed: correting MB images (2D phase, MUSE) at slice #' num2str(z) '/' num2str(Nz_MB) ', at tpn: ' num2str(tpn) '/' num2str(Nt_MB)]);
    end
end
