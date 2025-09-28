function [recon_image,PhaseMap_save,Phase_coff,xstart_end] = cpr2D_coils_mb2_SD_0930_2s(kspace,Image_SB_ref,CSP_SB,Phase2D_coeff,PhaseMap_1D,FOVshift,nshot)

sz = size(kspace); 
Ny_MB = sz(1);
Nx_MB = sz(2); 
Ncoil_MB = sz(3);

ns4 = nshot*2;

Image_ref_s0 = squeeze(Image_SB_ref(:,:,1,:));
Image_ref_s1 = squeeze(Image_SB_ref(:,:,2,:));

Image_pos_k1 = zeros(size(kspace));
Image_pos_k2 = zeros(size(kspace));
Image_neg_k3 = zeros(size(kspace));
Image_neg_k4 = zeros(size(kspace));

for c = 1 : Ncoil_MB
    kspace1 = zeros(Ny_MB,Nx_MB);
    kspace2 = zeros(Ny_MB,Nx_MB);
    kspace3 = zeros(Ny_MB,Nx_MB);
    kspace4 = zeros(Ny_MB,Nx_MB);
    
    kspace1(1:ns4:end,:) = kspace(1:ns4:end,:,c);
    kspace2(2:ns4:end,:) = kspace(2:ns4:end,:,c);
    kspace3(3:ns4:end,:) = kspace(3:ns4:end,:,c);
    kspace4(4:ns4:end,:) = kspace(4:ns4:end,:,c);
    
    Image_pos_k1(:,:,c) = fftshift(fft2(fftshift(kspace1)));
    Image_pos_k2(:,:,c) = fftshift(fft2(fftshift(kspace2)));
    Image_neg_k3(:,:,c) = fftshift(fft2(fftshift(kspace3)));
    Image_neg_k4(:,:,c) = fftshift(fft2(fftshift(kspace4)));
end

%% Create mask outside brain (background region)
RMS_tmp0 = RootMeanSquare(Image_ref_s0,1,Ncoil_MB,'Image');
[Mask0_wholebrain,~] = Mask_Generate(RMS_tmp0,0.16,0.16);

RMS_tmp1 = RootMeanSquare(Image_ref_s1,1,Ncoil_MB,'Image');
[Mask1_wholebrain,~] = Mask_Generate(RMS_tmp1,0.16,0.16);

Mask0_wholebrain_b = ~Mask0_wholebrain;
Mask1_wholebrain_b = ~Mask1_wholebrain;

xstart_s0 = find(sum(Mask0_wholebrain_b,1)>5,1);
xend_s0 = find(sum(Mask0_wholebrain_b,1)>5,1,'last');

xstart_s1 = find(sum(Mask1_wholebrain_b,1)>5,1);
xend_s1 = find(sum(Mask1_wholebrain_b,1)>5,1,'last');

xstart = min([xstart_s0 xstart_s1]);
xend = max([xend_s0 xend_s1]);

%% Correction along x direction
FOV = Ny_MB;

PhaseMap_s0 = zeros(Ny_MB,Nx_MB);
PhaseMap_s1 = zeros(Ny_MB,Nx_MB);

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

Phase_coff = zeros(Nx_MB,2,acc_MB);
Step_p1 = 0.02;
Step_p2 = 0.002;

for k1 = xstart : xend
    
    ref1_s0 = Phase2D_coeff(k1,1,3);
    ref2_s0 = Phase2D_coeff(k1,2,3);
    ref1_s1 = Phase2D_coeff(k1,1,3);
    ref2_s1 = Phase2D_coeff(k1,2,3);
    
    coeff1list_s0 = [ ref1_s0-Step_p1 ref1_s0 ref1_s0+Step_p1 ];
    coeff1ln_s0 = length(coeff1list_s0);
    coeff2list_s0 = [ ref2_s0-Step_p2 ref2_s0 ref2_s0+Step_p2 ];
    coeff2ln_s0 = length(coeff2list_s0);
    
    coeff1list_s1 = [ ref1_s1-Step_p1 ref1_s1 ref1_s1+Step_p1 ];
    coeff1ln_s1 = length(coeff1list_s1);
    coeff2list_s1 = [ ref2_s1-Step_p2 ref2_s1 ref2_s1+Step_p2 ];
    coeff2ln_s1 = length(coeff2list_s1);
     
    recon_line = zeros(coeff1ln_s0,coeff2ln_s0,coeff1ln_s1,coeff2ln_s1);
    
    for C1_s0 = 1 : coeff1ln_s0
        for C2_s0 = 1 : coeff2ln_s0
            for C1_s1 = 1 : coeff1ln_s1
                for C2_s1 = 1 : coeff2ln_s1
                            pc1_s0 = zeros(Nx_MB,1);
                            pc1_s0(:) = exp(1i*coeff1list_s0(C1_s0));
                            pc2_s0 = exp(1i*coeff2list_s0(C2_s0).*transpose(-Nx_MB/2:Nx_MB/2-1));
                            pc_s0 = pc1_s0.*pc2_s0.*exp(1i*(PhaseMap_1D(:,k1)));              
                            
                            pc1_s1 = zeros(Nx_MB,1);
                            pc1_s1(:) = exp(1i*coeff1list_s1(C1_s1));
                            pc2_s1 = exp(1i*coeff2list_s1(C2_s1).*transpose(-Nx_MB/2:Nx_MB/2-1));
                            pc_s1 = pc1_s1.*pc2_s1.*exp(1i*(PhaseMap_1D(:,k1)));                          
                            
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
                                
                                
                                pc_r = transpose([pc_s0(p_s0);pc_s1(p_s1)]);
                                
                                for c = 1 : Ncoil_MB
                                    CSP = transpose([CSP_SB(p_s0,k1,1,c);CSP_SB(p_s1,k1,2,c)]);
                                    
                                    CSP_matrix((k2-1)*Ncoil_MB*ns4 + (c-1)*ns4+1,(k2-1)*ns4*2 + 1:(k2-1)*ns4*2 + 8) = 1/ns4.*Pos_k1.*pc_r.*CSP;
                                    CSP_matrix((k2-1)*Ncoil_MB*ns4 + (c-1)*ns4+2,(k2-1)*ns4*2 + 1:(k2-1)*ns4*2 + 8) = 1/ns4.*Pos_k2.*pc_r.*CSP;
                                    
                                    CSP_matrix((k2-1)*Ncoil_MB*ns4 + (c-1)*ns4+3,(k2-1)*ns4*2 + 1:(k2-1)*ns4*2 + 8) = 1/ns4.*Neg_k3.*CSP;
                                    CSP_matrix((k2-1)*Ncoil_MB*ns4 + (c-1)*ns4+4,(k2-1)*ns4*2 + 1:(k2-1)*ns4*2 + 8) = 1/ns4.*Neg_k4.*CSP;
                                    
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
                                
                                LineCorr_s0(p_s0,1) = S((k2-1)*ns4*2 + 1:(k2-1)*ns4*2 + 4,1);
                                LineCorr_s1(p_s1,1) = S((k2-1)*ns4*2 + 5:(k2-1)*ns4*2 + 8,1);
                            end
                            
                            recon_line_tmp = abs(LineCorr_s0).*Mask0_wholebrain_b(:,k1) + abs(LineCorr_s1).*Mask1_wholebrain_b(:,k1);
                            recon_line(C1_s0,C2_s0,C1_s1,C2_s1) = squeeze(sum(recon_line_tmp,1));
                            
                end
            end
        end
    end
        
    tmp = find(recon_line == min(recon_line(:)));
    sz = [coeff1ln_s0 coeff2ln_s0 coeff1ln_s1 coeff2ln_s1];
    [p1_s0,p2_s0,p1_s1,p2_s1] = ind2sub(sz,tmp);
    
    Phase_coff(k1,:,1) = [coeff1list_s0(p1_s0)  coeff2list_s0(p2_s0)];
    Phase_coff(k1,:,2) = [coeff1list_s1(p1_s1)  coeff2list_s1(p2_s1)];    
    
    PhaseMap_s0(:,k1) = coeff1list_s0(p1_s0) + coeff2list_s0(p2_s0).*transpose(-Ny_MB/2:Ny_MB/2-1);
    PhaseMap_s1(:,k1) = coeff1list_s1(p1_s1) + coeff2list_s1(p2_s1).*transpose(-Ny_MB/2:Ny_MB/2-1);
    
    % parfor_progress;
    
end

%% Reconstruction
PhaseMap = zeros(Ny_MB,Nx_MB,acc_MB);
PhaseMap_SD = zeros(size(PhaseMap));
PhaseMap(:,:,1) = PhaseMap_s0;
PhaseMap(:,:,2) = PhaseMap_s1;

xstart_fit = [xstart_s0 xstart_s1];
xend_fit = [xend_s0 xend_s1];
xstart_end = [xstart_fit;xend_fit];
for c = 1 : acc_MB
    PhaseMap_unwrap2d = phaseunwrap2d(exp(1i.*(PhaseMap(:,:,c))));
    
    fmap2Dpre = PhaseMap_unwrap2d;
    
    xp = repmat((1:Ny_MB)',length(xstart_fit(c):xend_fit(c)),1);
    yp = reshape(repmat(xstart_fit(c):xend_fit(c),Ny_MB,1),[],1);
    
    PhaseMapS = fmap2Dpre(:,xstart_fit(c):xend_fit(c));
    zz = reshape(PhaseMapS,[],1);
    FitPara = fit([xp yp],zz,'poly23');
    
    xp_1 = repmat((1:Ny_MB)',Nx_MB,1);
    yp_1 = reshape(repmat(1:Nx_MB,Ny_MB,1),[],1);
    
    f = FitPara.p00 + FitPara.p10*xp_1 + FitPara.p01*yp_1 + FitPara.p20*xp_1.^2 + ...
        + FitPara.p11.*xp_1.*yp_1 + FitPara.p02.*yp_1.^2 + FitPara.p21.*xp_1.^2.*yp_1 + ...
        + FitPara.p12.*xp_1.*yp_1.^2 + FitPara.p03.*yp_1.^3;
    fmap2Dfit = reshape(f,Ny_MB,Nx_MB);
    
    PhaseMap_SD(:,:,c) = phaseunwrap2d(exp(1i.*(fmap2Dfit + PhaseMap_1D)));
end
PhaseMap_save = PhaseMap_SD;

recon_image = zeros(Ny_MB,Nx_MB,acc_MB);

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
        
        pc_r = exp(1i.*transpose([PhaseMap_SD(p_s0,k1,1);PhaseMap_SD(p_s1,k1,2)]));
        for c = 1 : Ncoil_MB
            CSP = transpose([CSP_SB(p_s0,k1,1,c);CSP_SB(p_s1,k1,2,c)]);
            CSP_matrix((k2-1)*Ncoil_MB*ns4 + (c-1)*ns4+1,(k2-1)*ns4*2 + 1:(k2-1)*ns4*2 + 8) = 1/ns4.*Pos_k1.*pc_r.*CSP;
            CSP_matrix((k2-1)*Ncoil_MB*ns4 + (c-1)*ns4+2,(k2-1)*ns4*2 + 1:(k2-1)*ns4*2 + 8) = 1/ns4.*Pos_k2.*pc_r.*CSP;
            
            CSP_matrix((k2-1)*Ncoil_MB*ns4 + (c-1)*ns4+3,(k2-1)*ns4*2 + 1:(k2-1)*ns4*2 + 8) = 1/ns4.*Neg_k3.*CSP;
            CSP_matrix((k2-1)*Ncoil_MB*ns4 + (c-1)*ns4+4,(k2-1)*ns4*2 + 1:(k2-1)*ns4*2 + 8) = 1/ns4.*Neg_k4.*CSP;
            
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
        
        LineCorr_s0(p_s0,1) = S((k2-1)*ns4*2 + 1:(k2-1)*ns4*2 + 4,1);
        LineCorr_s1(p_s1,1) = S((k2-1)*ns4*2 + 5:(k2-1)*ns4*2 + 8,1);
    end
    
    recon_image(:,k1,1) = LineCorr_s0;
    recon_image(:,k1,2) = LineCorr_s1;
    
end

end



