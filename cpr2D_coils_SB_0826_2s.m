function [Phase_coeff,CorrImg,PhaseMap_pc] = cpr2D_coils_SB_0826_2s(kspace,Image_ref,CSP_SB,PhaseMap_1D,Nshot,weight_csp)

sz = size(kspace); 
Ny_SB = sz(1);
Nx_SB = sz(2);
Ncoil_SB = sz(3);

ns4 = Nshot*2;

Image_ref_s0 = Image_ref;

Image_pos_k1 = zeros(size(kspace));
Image_pos_k2 = zeros(size(kspace));
Image_neg_k3 = zeros(size(kspace));
Image_neg_k4 = zeros(size(kspace));

for c = 1 : Ncoil_SB
    kspace1 = zeros(Ny_SB,Nx_SB);
    kspace2 = zeros(Ny_SB,Nx_SB);
    kspace3 = zeros(Ny_SB,Nx_SB);
    kspace4 = zeros(Ny_SB,Nx_SB);
    
    kspace1(1:ns4:end,:) = kspace(1:ns4:end,:,c);
    kspace2(2:ns4:end,:) = kspace(2:ns4:end,:,c);
    kspace3(3:ns4:end,:) = kspace(3:ns4:end,:,c);
    kspace4(4:ns4:end,:) = kspace(4:ns4:end,:,c);
    
    Image_pos_k1(:,:,c) = fftshift(fft2(fftshift(kspace1)));
    Image_pos_k2(:,:,c) = fftshift(fft2(fftshift(kspace2)));
    Image_neg_k3(:,:,c) = fftshift(fft2(fftshift(kspace3)));
    Image_neg_k4(:,:,c) = fftshift(fft2(fftshift(kspace4)));
end

%% Create mask
RMS_tmp0 = RootMeanSquare(Image_ref_s0,1,Ncoil_SB,'Image');
[Mask_s0,~] = Mask_Generate(RMS_tmp0,0.16,0.16);

tmp = Mask_s0;
tmp1 = tmp > 0 ;
a1 = tmp1 + Mask_s0([Nx_SB/2+1:end 1:Nx_SB/2],:);
a1 = a1 >1;
Mask_csp = (a1 > 0);

xstart_csp = find(sum(Mask_csp,1)>0,1);
xend_csp = find(sum(Mask_csp,1)>0,1,'last');

Mask_bg = zeros(Ny_SB,Nx_SB);
tmp2 = tmp > 0;
xstart_bg = find(sum(tmp2,1)>5,1);
xend_bg = find(sum(tmp2,1)>5,1,'last');
Mask_bg(:,sum(tmp2,1)>0) = 1;
Mask_bg = Mask_bg - tmp2;
Mask_bg(:,[1:xstart_bg-1 xend_bg+1:end]) = 0;
mskld_bg = sum(Mask_bg,1);

xstart = min([xstart_bg xstart_csp]);
xend = max([xend_bg xend_csp]);

%% Determine central column
tmp = sum(Mask_s0,1);
CenColumn_p = find(tmp == max(tmp(:)));
if length(CenColumn_p) > 1
    CenColumn_p = CenColumn_p(ceil(length(CenColumn_p)/2));
end
xlocation0 = CenColumn_p;

%% Correction along x direction
Pos_k1 = [1 1 1 1];
Pos_k2 = [1 exp(1i*pi/2) exp(1i*pi) exp(-1i*pi/2)];
Neg_k3 = [1 -1 1 -1];
Neg_k4 = [1 exp(-1i*pi/2) exp(-1i*pi) exp(1i*pi/2)];

FOV = Ny_SB;
Phase_coeff = zeros(Nx_SB,2);
coeff1step_c2 = 0.1;
coeff2step_c2 = 0.01;
savefactor1 = 0.15;
savefactor2 = 0.015;

Difference = zeros(Nx_SB,2,2);
ImageCorr = zeros(Ny_SB,Nx_SB);

num_search = 1;
Index_csp = zeros(Nx_SB,2,2);
for k1 = xlocation0 : xend
    if k1 == xlocation0
        coeff1list = round(-savefactor1*pi,1) : coeff1step_c2 : round(savefactor1*pi,1);
        coeff1ln = length(coeff1list);
        coeff2list = round(-savefactor2*pi,2) : coeff2step_c2 : round(savefactor2*pi,2);
        coeff2ln = length(coeff2list);
    else
        coeff1list = ref1 - num_search*coeff1step_c2 : coeff1step_c2 : ref1 + num_search*coeff1step_c2;
        coeff1ln = length(coeff1list);
        coeff2list = ref2 - num_search*coeff2step_c2 : coeff2step_c2 : ref2 + num_search*coeff2step_c2;
        coeff2ln = length(coeff2list);
    end
    
    LineCorr = zeros(Ny_SB,coeff1ln,coeff2ln,Ncoil_SB);
    DiffSum_csp = zeros(coeff1ln,coeff2ln,Ncoil_SB);
    DiffSum_bg = zeros(coeff1ln,coeff2ln);
    LineCorrCombine = zeros(Ny_SB,coeff1ln,coeff2ln);
    CoilSensProf_Ite = zeros(Ny_SB,coeff1ln,coeff2ln,Ncoil_SB);
    
    for C1 =  1 : coeff1ln  % x dirction
        for C2 = 1 : coeff2ln % y dirction
            
            pc1 = coeff1list(C1);
            pc2 = transpose(-Ny_SB/2 : Ny_SB/2-1);
            pc2 = coeff2list(C2)*pc2;
            pc = pc1 + pc2 + PhaseMap_1D(:,k1) ;
            
            CSP_matrix =  zeros(ns4*Ncoil_SB*Ny_SB/ns4,Ny_SB*Ncoil_SB); % csp and phase matrix
            Imat = zeros(ns4*Ncoil_SB*Ny_SB/ns4,1); % Odd-even matrix
            
            for k2 = 1 : Ny_SB/ns4  % y direction
                p = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4];
                p(p>FOV) = p(p>FOV) - FOV;
                
                pc_r = exp(1i.*transpose(pc(p,1)));
                
                for c = 1 : Ncoil_SB
                    CSP_matrix((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4+1,(k2-1)*Ncoil_SB*ns4 + (c-1)*4+1:(k2-1)*Ncoil_SB*ns4 + 4*c) = 1/ns4.*Pos_k1.*pc_r;
                    CSP_matrix((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4+2,(k2-1)*Ncoil_SB*ns4 + (c-1)*4+1:(k2-1)*Ncoil_SB*ns4 + 4*c) = 1/ns4.*Pos_k2.*pc_r;
                    
                    CSP_matrix((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4+3,(k2-1)*Ncoil_SB*ns4 + (c-1)*4+1:(k2-1)*Ncoil_SB*ns4 + 4*c) = 1/ns4.*Neg_k3;
                    CSP_matrix((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4+4,(k2-1)*Ncoil_SB*ns4 + (c-1)*4+1:(k2-1)*Ncoil_SB*ns4 + 4*c) = 1/ns4.*Neg_k4;
                    
                end
                Imat((k2-1)*Ncoil_SB*ns4 + 1 : ns4 : k2*Ncoil_SB*ns4) = squeeze(Image_pos_k1(k2,k1,:));
                Imat((k2-1)*Ncoil_SB*ns4 + 2 : ns4 : k2*Ncoil_SB*ns4) = squeeze(Image_pos_k2(k2,k1,:));
                Imat((k2-1)*Ncoil_SB*ns4 + 3 : ns4 : k2*Ncoil_SB*ns4) = squeeze(Image_neg_k3(k2,k1,:));
                Imat((k2-1)*Ncoil_SB*ns4 + 4 : ns4 : k2*Ncoil_SB*ns4) = squeeze(Image_neg_k4(k2,k1,:));
            end
            
            Phase = CSP_matrix;
            S = Phase\Imat;
            
            for k2 = 1 : Ny_SB/ns4
                p = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4];
                p(p>FOV) = p(p>FOV) - FOV;
                for c = 1 : Ncoil_SB
                    LineCorr(p,C1,C2,c) = S((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4 + 1 : (k2-1)*Ncoil_SB*ns4 + (c-1)*ns4 + 4);
                end
            end
            
            LineCorrCombine(:,C1,C2) = squeeze(sqrt(mean(abs(LineCorr(:,C1,C2,:)).^2,4))); % Combine coils
            
            CoilSensProf_Ite(:,C1,C2,:) = squeeze(LineCorr(:,C1,C2,:)./LineCorrCombine(:,C1,C2));
            
            CoilSensProf_Ite(:,C1,C2,:) = squeeze(CoilSensProf_Ite(:,C1,C2,:)./CoilSensProf_Ite(:,C1,C2,1).*abs(CoilSensProf_Ite(:,C1,C2,1)));
            
            for c = 1 : Ncoil_SB
                q1_s1 = CoilSensProf_Ite(:,C1,C2,c).*Mask_csp(:,k1);
                q2_s1 = CSP_SB(:,k1,c).*Mask_csp(:,k1);
                Diff_s1 = abs(abs(q1_s1) - abs(q2_s1));
                DiffSum_csp(C1,C2,c) = sum(Diff_s1,1);                
            end
            if mskld_bg(k1) ~=0
                DiffSum_bg(C1,C2) = squeeze(sum(LineCorrCombine(:,C1,C2).*Mask_bg(:,k1),1));
            end            
        end
    end
    
    if sum(Mask_csp(:,k1)) ~= 0
        tmp1 = squeeze(sqrt(mean(abs(DiffSum_csp).^2,3))); % Combine coils;
        [p1_csp,p2_csp] = find(tmp1 == min(tmp1(:)));
        Index_csp(k1,1,:) = [p1_csp,p2_csp];
        if mskld_bg(k1) ~= 0 % csp ~= 0, bg~= 0
            [p1_bg,p2_bg] = find(DiffSum_bg == min(DiffSum_bg(:)));
            ImageCorr(:,k1) = LineCorrCombine(:,p1_bg,p2_bg);
            Index_csp(k1,2,:) = [p1_bg,p2_bg];
            Phase_coeff(k1,:,1) = [coeff1list(p1_csp) coeff2list(p2_csp)];
            Phase_coeff(k1,:,2) = [coeff1list(p1_bg) coeff2list(p2_bg)];
            Phase_coeff(k1,:,3) = [coeff1list(p1_csp).*weight_csp + coeff1list(p1_bg).*(1 - weight_csp) coeff2list(p2_csp).*weight_csp + coeff2list(p2_bg).*(1 - weight_csp)];
            Difference(k1,1,:) = [min(tmp1(:)) max(tmp1(:))];
            Difference(k1,2,:) = [min(DiffSum_bg(:)) max(DiffSum_bg(:))];
        else % csp ~= 0, bg = 0
            Phase_coeff(k1,:,1) = [coeff1list(p1_csp) coeff2list(p2_csp)];
            Phase_coeff(k1,:,2) = [coeff1list(p1_csp) coeff2list(p2_csp)];
            Phase_coeff(k1,:,3) = [coeff1list(p1_csp) coeff2list(p2_csp)];
            Difference(k1,1,:) = [min(tmp1(:)) max(tmp1(:))];
            Difference(k1,2,:) = [0 0];
            
        end
    else
        if mskld_bg(k1) ~= 0 % csp = 0, bg~= 0
            [p1_bg,p2_bg] = find(DiffSum_bg == min(DiffSum_bg(:)));
            ImageCorr(:,k1) = LineCorrCombine(:,p1_bg,p2_bg);
            Index_csp(k1,2,:) = [p1_bg,p2_bg];
            Phase_coeff(k1,:,1) = [coeff1list(p1_bg) coeff2list(p2_bg)];
            Phase_coeff(k1,:,2) = [coeff1list(p1_bg) coeff2list(p2_bg)];
            Phase_coeff(k1,:,3) = [coeff1list(p1_bg).*weight_csp + coeff1list(p1_bg).*(1 - weight_csp) coeff2list(p2_bg).*weight_csp + coeff2list(p2_bg).*(1 - weight_csp)];
            
            Difference(k1,1,:) = [0 0];
            Difference(k1,2,:) = [min(DiffSum_bg(:)) max(DiffSum_bg(:))];
        end        
    end
    
    ref1 = Phase_coeff(k1,1,3);
    ref2 = Phase_coeff(k1,2,3);

    if k1 == xlocation0
        ref001 = ref1;
        ref002 = ref2;
    end
    
    % parfor_progress;
    
end
% parfor_progress(0);

%
ref1 = ref001;
ref2 = ref002;
% parfor_progress(length(xlocation0-1:-1:xstart));
for k1 = xlocation0-1 :-1:xstart
    
    coeff1list = ref1 - num_search*coeff1step_c2 : coeff1step_c2 : ref1 + num_search*coeff1step_c2;
    coeff1ln = length(coeff1list);
    coeff2list = ref2 - num_search*coeff2step_c2 : coeff2step_c2 : ref2 + num_search*coeff2step_c2;
    coeff2ln = length(coeff2list);     
    
    LineCorr = zeros(Ny_SB,coeff1ln,coeff2ln,Ncoil_SB);
    DiffSum_csp = zeros(coeff1ln,coeff2ln,Ncoil_SB);
    DiffSum_bg = zeros(coeff1ln,coeff2ln);
    LineCorrCombine = zeros(Ny_SB,coeff1ln,coeff2ln);
    CoilSensProf_Ite = zeros(Ny_SB,coeff1ln,coeff2ln,Ncoil_SB);
    
    for C1 =  1 : coeff1ln  % x dirction
        for C2 = 1 : coeff2ln % y dirction
            
            pc1 = coeff1list(C1);
            pc2 = transpose(-Ny_SB/2 : Ny_SB/2-1);
            pc2 = coeff2list(C2)*pc2;
            pc = pc1 + pc2 + PhaseMap_1D(:,k1) ;
            
            CSP_matrix =  zeros(ns4*Ncoil_SB*Ny_SB/ns4,Ny_SB*Ncoil_SB); % csp and phase matrix
            Imat = zeros(ns4*Ncoil_SB*Ny_SB/ns4,1); % Odd-even matrix
            
            for k2 = 1 : Ny_SB/ns4  % y direction
                p = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4];
                p(p>FOV) = p(p>FOV) - FOV;
                
                pc_r = exp(1i.*transpose(pc(p,1)));
                
                for c = 1 : Ncoil_SB
                    CSP_matrix((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4+1,(k2-1)*Ncoil_SB*ns4 + (c-1)*4+1:(k2-1)*Ncoil_SB*ns4 + 4*c) = 1/ns4.*Pos_k1.*pc_r;
                    CSP_matrix((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4+2,(k2-1)*Ncoil_SB*ns4 + (c-1)*4+1:(k2-1)*Ncoil_SB*ns4 + 4*c) = 1/ns4.*Pos_k2.*pc_r;
                    
                    CSP_matrix((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4+3,(k2-1)*Ncoil_SB*ns4 + (c-1)*4+1:(k2-1)*Ncoil_SB*ns4 + 4*c) = 1/ns4.*Neg_k3;
                    CSP_matrix((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4+4,(k2-1)*Ncoil_SB*ns4 + (c-1)*4+1:(k2-1)*Ncoil_SB*ns4 + 4*c) = 1/ns4.*Neg_k4;
                    
                end
                Imat((k2-1)*Ncoil_SB*ns4 + 1 : ns4 : k2*Ncoil_SB*ns4) = squeeze(Image_pos_k1(k2,k1,:));
                Imat((k2-1)*Ncoil_SB*ns4 + 2 : ns4 : k2*Ncoil_SB*ns4) = squeeze(Image_pos_k2(k2,k1,:));
                Imat((k2-1)*Ncoil_SB*ns4 + 3 : ns4 : k2*Ncoil_SB*ns4) = squeeze(Image_neg_k3(k2,k1,:));
                Imat((k2-1)*Ncoil_SB*ns4 + 4 : ns4 : k2*Ncoil_SB*ns4) = squeeze(Image_neg_k4(k2,k1,:));
            end
            
            Phase = CSP_matrix;
            S = Phase\Imat;
            
            for k2 = 1 : Ny_SB/ns4
                p = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4];
                p(p>FOV) = p(p>FOV) - FOV;
                for c = 1 : Ncoil_SB
                    LineCorr(p,C1,C2,c) = S((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4 + 1 : (k2-1)*Ncoil_SB*ns4 + (c-1)*ns4 + 4);
                end
            end
            
            LineCorrCombine(:,C1,C2) = squeeze(sqrt(mean(abs(LineCorr(:,C1,C2,:)).^2,4))); % Combine coils
            
            CoilSensProf_Ite(:,C1,C2,:) = squeeze(LineCorr(:,C1,C2,:)./LineCorrCombine(:,C1,C2));
            
            CoilSensProf_Ite(:,C1,C2,:) = squeeze(CoilSensProf_Ite(:,C1,C2,:)./CoilSensProf_Ite(:,C1,C2,1).*abs(CoilSensProf_Ite(:,C1,C2,1)));
            
            for c = 1 : Ncoil_SB
                q1_s1 = CoilSensProf_Ite(:,C1,C2,c).*Mask_csp(:,k1);
                q2_s1 = CSP_SB(:,k1,c).*Mask_csp(:,k1);
                Diff_s1 = abs(abs(q1_s1) - abs(q2_s1));
                DiffSum_csp(C1,C2,c) = sum(Diff_s1,1);
                
            end
            if mskld_bg(k1) ~=0
                DiffSum_bg(C1,C2) = squeeze(sum(LineCorrCombine(:,C1,C2).*Mask_bg(:,k1),1));
            end            
        end
    end
    
    if sum(Mask_csp(:,k1)) ~= 0
        tmp1 = squeeze(sqrt(mean(abs(DiffSum_csp).^2,3))); % Combine coils;
        [p1_csp,p2_csp] = find(tmp1 == min(tmp1(:)));
        Index_csp(k1,1,:) = [p1_csp,p2_csp];
        if mskld_bg(k1) ~= 0 % csp ~= 0, bg~= 0
            [p1_bg,p2_bg] = find(DiffSum_bg == min(DiffSum_bg(:)));
            ImageCorr(:,k1) = LineCorrCombine(:,p1_bg,p2_bg);
            Index_csp(k1,2,:) = [p1_bg,p2_bg];
            Phase_coeff(k1,:,1) = [coeff1list(p1_csp) coeff2list(p2_csp)];
            Phase_coeff(k1,:,2) = [coeff1list(p1_bg) coeff2list(p2_bg)];
            Phase_coeff(k1,:,3) = [coeff1list(p1_csp).*weight_csp + coeff1list(p1_bg).*(1 - weight_csp) coeff2list(p2_csp).*weight_csp + coeff2list(p2_bg).*(1 - weight_csp)];
            Difference(k1,1,:) = [min(tmp1(:)) max(tmp1(:))];
            Difference(k1,2,:) = [min(DiffSum_bg(:)) max(DiffSum_bg(:))];
        else % csp ~= 0, bg = 0
            Phase_coeff(k1,:,1) = [coeff1list(p1_csp) coeff2list(p2_csp)];
            Phase_coeff(k1,:,2) = [coeff1list(p1_csp) coeff2list(p2_csp)];
            Phase_coeff(k1,:,3) = [coeff1list(p1_csp) coeff2list(p2_csp)];
            Difference(k1,1,:) = [min(tmp1(:)) max(tmp1(:))];
            Difference(k1,2,:) = [0 0];
            
        end
    else
        if mskld_bg(k1) ~= 0 % csp = 0, bg~= 0
            [p1_bg,p2_bg] = find(DiffSum_bg == min(DiffSum_bg(:)));
            ImageCorr(:,k1) = LineCorrCombine(:,p1_bg,p2_bg);
            Index_csp(k1,2,:) = [p1_bg,p2_bg];
            Phase_coeff(k1,:,1) = [coeff1list(p1_bg) coeff2list(p2_bg)];
            Phase_coeff(k1,:,2) = [coeff1list(p1_bg) coeff2list(p2_bg)];
            Phase_coeff(k1,:,3) = [coeff1list(p1_bg).*weight_csp + coeff1list(p1_bg).*(1 - weight_csp) coeff2list(p2_bg).*weight_csp + coeff2list(p2_bg).*(1 - weight_csp)];
            
            Difference(k1,1,:) = [0 0];
            Difference(k1,2,:) = [min(DiffSum_bg(:)) max(DiffSum_bg(:))];
        end        
    end   

    ref1 = Phase_coeff(k1,1,3);
    ref2 = Phase_coeff(k1,2,3);
    
    % parfor_progress;
    
end
% parfor_progress(0);

%% Full matrix inversion for correction
acc_MB = 1;

PhaseMap = zeros(Ny_SB,Nx_SB);
for k1 = xstart : xend
    PhaseMap(:,k1) = exp(1i.*Phase_coeff(k1,2,3).*transpose(-Nx_SB/2:Nx_SB/2-1)).*exp(1i*repmat(Phase_coeff(k1,1,3),Ny_SB,1));
end

PhaseMap_unwrap2d = phaseunwrap2d(PhaseMap);

fmap2Dpre = PhaseMap_unwrap2d;

xp = repmat((1:Ny_SB)',length(xstart:xend),1);
yp = reshape(repmat(xstart:xend,Ny_SB,1),[],1);

PhaseMapS = fmap2Dpre(:,xstart:xend);
zz = reshape(PhaseMapS,[],1);
FitPara = fit([xp yp],zz,'poly23');

xp_1 = repmat((1:Ny_SB)',Nx_SB,1);
yp_1 = reshape(repmat(1:Nx_SB,Ny_SB,1),[],1);

f = FitPara.p00 + FitPara.p10*xp_1 + FitPara.p01*yp_1 + FitPara.p20*xp_1.^2 + ...
    + FitPara.p11.*xp_1.*yp_1 + FitPara.p02.*yp_1.^2 + FitPara.p21.*xp_1.^2.*yp_1 + ...
    + FitPara.p12.*xp_1.*yp_1.^2 + FitPara.p03.*yp_1.^3;
fmap2Dfit = reshape(f,Ny_SB,Nx_SB);

PhaseMap_pc = phaseunwrap2d(exp(1i.*(fmap2Dfit + PhaseMap_1D)));

recon_image = zeros(Ny_SB,Nx_SB,acc_MB);

for k1 = 1 : Nx_SB % x direction
    LineCorr_s1 = zeros(Ny_SB,1);
    
    CSP_matrix =  zeros(ns4*Ncoil_SB*Ny_SB/ns4,Ny_SB); % csp and phase matrix
    Imat = zeros(ns4*Ncoil_SB*Ny_SB/ns4,1); % Odd-even matrix
    
    for k2 = 1 : Ny_SB/ns4  % y direction
        p = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4];
        p(p>FOV) = p(p>FOV) - FOV;
        
        pc_r = exp(1i.*transpose(PhaseMap_pc(p,k1)));
        
        for c = 1 : Ncoil_SB
            CSP = transpose(CSP_SB(p,k1,c));
            CSP_matrix((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4+1,(k2-1)*ns4 + 1:(k2-1)*ns4 + 4) = 1/ns4.*Pos_k1.*pc_r.*CSP;
            CSP_matrix((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4+2,(k2-1)*ns4 + 1:(k2-1)*ns4 + 4) = 1/ns4.*Pos_k2.*pc_r.*CSP;
            
            CSP_matrix((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4+3,(k2-1)*ns4 + 1:(k2-1)*ns4 + 4) = 1/ns4.*Neg_k3.*CSP;
            CSP_matrix((k2-1)*Ncoil_SB*ns4 + (c-1)*ns4+4,(k2-1)*ns4 + 1:(k2-1)*ns4 + 4) = 1/ns4.*Neg_k4.*CSP;
            
        end
        Imat((k2-1)*Ncoil_SB*ns4 + 1 : ns4 : k2*Ncoil_SB*ns4) = squeeze(Image_pos_k1(k2,k1,:));
        Imat((k2-1)*Ncoil_SB*ns4 + 2 : ns4 : k2*Ncoil_SB*ns4) = squeeze(Image_pos_k2(k2,k1,:));
        Imat((k2-1)*Ncoil_SB*ns4 + 3 : ns4 : k2*Ncoil_SB*ns4) = squeeze(Image_neg_k3(k2,k1,:));
        Imat((k2-1)*Ncoil_SB*ns4 + 4 : ns4 : k2*Ncoil_SB*ns4) = squeeze(Image_neg_k4(k2,k1,:));
    end
    
    Phase = CSP_matrix;
    S = Phase\Imat;
    
    for k2 = 1 : Ny_SB/ns4
        p = [k2;k2+FOV/ns4;k2+FOV/2;k2+FOV*3/4];
        p(p>FOV) = p(p>FOV) - FOV;
        
        LineCorr_s1(p,1) = S((k2-1)*ns4 + 1:(k2-1)*ns4 + 4,1);

    end    
    
    recon_image(:,k1,1) = LineCorr_s1;
end

CorrImg = recon_image;


end

