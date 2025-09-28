% % Use average image to perform correction -> phase map ?2shot correction?
% Then apply to all slices and all tpn
kspace_MB_2shot = kdataoriginal;
MB_SliceNum = [1: Nz_MB; 1 + Nz_MB: Nz_MB*2];

%% Generate mask and determine slice
% Normalize
imgtmp0 = Image_SB_ref*3e4/max(RMS_SB_ref(:)); % To find a mask
RMS_EPI_pro = RootMeanSquare(imgtmp0,Nz_SB,Ncoil_MB,'Image');

% Select slice based on mask area
hanning20 = hanning(floor(Nx_MB/8))*hanning(floor(Nx_MB/8))';
quickrecon2 = RMS_EPI_pro*0;
for z = 1 : Nz_MB    
    tmp0 = RMS_EPI_pro(:,:,MB_SliceNum(1,z));
    tmp1 = RMS_EPI_pro(:,:,MB_SliceNum(2,z));  % Select one image

    quickrecon2(:,:,MB_SliceNum(1,z)) = conv2(tmp0,hanning20,'same'); % Filtered images   
    quickrecon2(:,:,MB_SliceNum(2,z)) = conv2(tmp1,hanning20,'same'); % Filtered images   
end
quickrecon2mask = quickrecon2 > 0.1*max(quickrecon2(:));  % Create a mask
clear tmp1 tmp2;
tmp2 = zeros(Nz_SB,1);

for z = 1: Nz_SB
    tmp1 = quickrecon2mask(:,:,z);
    tmp2(z) = length(find(tmp1)); % Find the number of 1 in mask
end
tmp2 = tmp2(MB_SliceNum(1,:)) + tmp2(MB_SliceNum(2,:));
[tmp3,tmp3i] = sort(tmp2,'descend');
tmp4 = find(tmp3>0.8*tmp3(1));
tmp5 = tmp4(1);
tmp6 = tmp3i(tmp5);
pick_a_slice = tmp6;
choosethisslice = pick_a_slice ;

%% Determine mask of different regions: csp and bg
RMS_tmp0 = RMS_SB_ref(:,:,MB_SliceNum(1,choosethisslice));
[Mask_s0,Mask_threshold_MB_s0] = Mask_Generate(RMS_tmp0,0.16,0.16);

RMS_tmp1 = RMS_SB_ref(:,:,MB_SliceNum(2,choosethisslice)); 
[Mask_s1,Mask_threshold_MB_s1] = Mask_Generate(RMS_tmp1,0.16,0.16);

% if FOVshift == 0
%     if nshot == 1
%         tmp = Mask_s0 + Mask_s1;
%         tmp1 = tmp > 0 ;
%         a1 = tmp1 + Mask_s0([Ny_MB/2+1:end 1:Ny_MB/2],:);
%         a1 = a1 >1;
%         a2 = tmp1 + Mask_s1([Ny_MB/2+1:end 1:Ny_MB/2],:);
%         a2 = a2 >1;
% 
%         Mask_csp = (a1 + a2 > 0);
%         Mask_csp = imfill(Mask_csp,'holes');
%         Mask_csp = bwareaopen(Mask_csp,ceil(Ny_MB*0.1));
%         xstart_csp = find(sum(Mask_csp,1)>5,1);
%         xend_csp = find(sum(Mask_csp,1)>5,1,'last');
% 
%         Mask_bg = zeros(Ny_MB,Nx_MB);
%         tmp2 = tmp > 0;
%         xstart_bg = find(sum(tmp2,1)>5,1);
%         xend_bg = find(sum(tmp2,1)>5,1,'last');
%         Mask_bg(:,sum(tmp2,1)>0) = 1;
%         Mask_bg = Mask_bg - tmp2;
%         Mask_bg(:,[1:xstart_bg-1 xend_bg+1:end]) = 0;
%         Mask_bg = bwareaopen(Mask_bg,ceil(Ny_MB*0.1));
%         mskld_bg = sum(Mask_bg,1);
% 
%         Mask_background = zeros(Ny_SB,Nx_SB);
%         Mask_background(:,[1:xstart_bg-1 xend_bg+1:end]) = 1;
%     end
% 
%     if nshot == 2
%         tmp = Mask_s0 + Mask_s1;
%         tmp1 = tmp > 0 ;
%         a1 = tmp1 + Mask_s0([Ny_SB/4+1:end 1:Ny_SB/4],:);
%         a1 = a1 > 1;
%         a2 = tmp1 + Mask_s0([Ny_SB*3/4+1:end 1:Ny_SB*3/4],:);
%         a2 = a2 >1;
% 
%         a3 = tmp1 + Mask_s1([Ny_SB/4+1:end 1:Ny_SB/4],:);
%         a3 = a3 > 1;
%         a4 = tmp1 + Mask_s1([Ny_SB*3/4+1:end 1:Ny_SB*3/4],:);
%         a4 = a4 >1;
%         Mask_csp = (a1 + a2 + a3 + a4) > 0;
% 
%         Mask_csp = imfill(Mask_csp,'holes');
%         Mask_csp = bwareaopen(Mask_csp,ceil(Ny_SB*0.1));
%         xstart_csp = find(sum(Mask_csp,1)>5,1);
%         xend_csp = find(sum(Mask_csp,1)>5,1,'last');
% 
%         Mask_bg = zeros(Ny_SB,Nx_SB);
%         tmp2 = tmp > 0;
%         xstart_bg = find(sum(tmp2,1)>5,1);
%         xend_bg = find(sum(tmp2,1)>5,1,'last');
%         Mask_bg(:,sum(tmp2,1)>0) = 1;
%         Mask_bg = Mask_bg - tmp2;
%         Mask_bg(:,[1:xstart_bg-1 xend_bg+1:end]) = 0;
%         Mask_bg = bwareaopen(Mask_bg,ceil(Ny_SB*0.1));
%         mskld_bg = sum(Mask_bg,1);
% 
%         Mask_background = zeros(Ny_SB,Nx_SB);
%         Mask_background(:,[1:xstart_bg-1 xend_bg+1:end]) = 1;
%     end
% end

if FOVshift == 1/2
   Mask_s1 = Mask_s1([Ny_MB/2+1:end 1:Ny_MB/2],:);% With FOV/2 shift

    if nshot == 1
        tmp = Mask_s0 + Mask_s1;
        tmp1 = tmp > 0 ;
        a1 = tmp1 + Mask_s0([Ny_MB/2+1:end 1:Ny_MB/2],:);
        a1 = a1 >1;
        a2 = tmp1 + Mask_s1([Ny_MB/2+1:end 1:Ny_MB/2],:);
        a2 = a2 >1;
        
        Mask_csp = (a1 + a2 > 0);
        Mask_csp = imfill(Mask_csp,'holes');
        Mask_csp = bwareaopen(Mask_csp,ceil(Ny_MB*0.1));
        xstart_csp = find(sum(Mask_csp,1)>5,1);
        xend_csp = find(sum(Mask_csp,1)>5,1,'last');
        
        Mask_bg = zeros(Ny_MB,Nx_MB);
        tmp2 = tmp > 0;
        xstart_bg = find(sum(tmp2,1)>5,1);
        xend_bg = find(sum(tmp2,1)>5,1,'last');
        Mask_bg(:,sum(tmp2,1)>0) = 1;
        Mask_bg = Mask_bg - tmp2;
        Mask_bg(:,[1:xstart_bg-1 xend_bg+1:end]) = 0;
        Mask_bg = bwareaopen(Mask_bg,ceil(Ny_MB*0.1));
        mskld_bg = sum(Mask_bg,1);
        
        Mask_background = zeros(Ny_SB,Nx_SB);
        Mask_background(:,[1:xstart_bg-1 xend_bg+1:end]) = 1;
    end
    
    if nshot == 2
        tmp = Mask_s0 + Mask_s1;
        tmp1 = tmp > 0 ;
        a1 = tmp1 + Mask_s0([Ny_SB/4+1:end 1:Ny_SB/4],:);
        a1 = a1 > 1;
        a2 = tmp1 + Mask_s0([Ny_SB*3/4+1:end 1:Ny_SB*3/4],:);
        a2 = a2 >1;
        
        a3 = tmp1 + Mask_s1([Ny_SB/4+1:end 1:Ny_SB/4],:);
        a3 = a3 > 1;
        a4 = tmp1 + Mask_s1([Ny_SB*3/4+1:end 1:Ny_SB*3/4],:);
        a4 = a4 >1;
        Mask_csp = (a1 + a2 + a3 + a4) > 0;
        
        Mask_csp = imfill(Mask_csp,'holes');
        Mask_csp = bwareaopen(Mask_csp,ceil(Ny_SB*0.1));
        xstart_csp = find(sum(Mask_csp,1)>5,1);
        xend_csp = find(sum(Mask_csp,1)>5,1,'last');
        
        Mask_bg = zeros(Ny_SB,Nx_SB);
        tmp2 = tmp > 0;
        xstart_bg = find(sum(tmp2,1)>5,1);
        xend_bg = find(sum(tmp2,1)>5,1,'last');
        Mask_bg(:,sum(tmp2,1)>0) = 1;
        Mask_bg = Mask_bg - tmp2;
        Mask_bg(:,[1:xstart_bg-1 xend_bg+1:end]) = 0;
        Mask_bg = bwareaopen(Mask_bg,ceil(Ny_SB*0.1));
        mskld_bg = sum(Mask_bg,1);
        
        Mask_background = zeros(Ny_SB,Nx_SB);
        Mask_background(:,[1:xstart_bg-1 xend_bg+1:end]) = 1;
    end
end

%% Run CPR on one slice: cpr1D linear
% Define cycled phase step and range
stepN1 = 31;
stepN2 = 31;
SaveFactor1 = 0.05;
SaveFactor2 = 0.1;
coeff1step = SaveFactor1*2*pi/(stepN1-1);
coeff1list = -pi*SaveFactor1 : coeff1step : pi*SaveFactor1;
coeff2step = SaveFactor2*2*pi/(stepN2-1);
coeff2list = transpose((-(SaveFactor2*pi):coeff2step:(SaveFactor2*pi)));
coeff2ln = length(coeff2list);
pc1 = exp(1i*coeff1list); % Intersect
pc2 = -Nx_SB/2:Nx_SB/2-1;
tmp = (2*coeff2list./Ny_SB)*pc2;
pc2 = exp(1i*tmp); % Slope

CSP_tmp_pc = zeros(Ny_MB,Nx_MB,length(pc1),size(pc2,1),Ncoil_MB);

IC_tmp_pc = zeros(Ny_MB,Nx_MB,length(pc1),size(pc2,1),Ncoil_MB);

for k1 = 1: length(pc1)
    for k2 = 1 : size(pc2,1)
        IC_cpr1D_oneslice = zeros(Ny_MB,Nx_MB,Ncoil_MB);
        kspace_cpr1D_oneslice =  zeros(Ny_MB,Nx_MB,Ncoil_MB);
        for c = 1 : Ncoil_MB
            for y = 1 : nshot                
                arrtmp = kspace_MB_2shot(y+nshot:ns4:end,:,choosethisslice,c)*pc1(k1);
                ktmp = fftshift(fft(fftshift(arrtmp,2),[],2),2);
                pc2tmp = pc2(k2,:);
                pc2tmp = repmat(pc2tmp,[Ny_MB/ns4 1]);
                kcorr = ktmp.*pc2tmp;
                kspace_cpr1D_oneslice(y+nshot:ns4:end,:,c) = ifftshift(ifft(ifftshift(kcorr,2),[],2),2);
                kspace_cpr1D_oneslice(y:ns4:end,:,c) = kspace_MB_2shot(y:ns4:end,:,choosethisslice,c);
            end
            IC_cpr1D_oneslice(:,:,c) = fftshift(fft2(fftshift(kspace_cpr1D_oneslice(:,:,c))));
        end
        
        RMStmp = RootMeanSquare(IC_cpr1D_oneslice,1,Ncoil_MB,'Image');
        CSP_tmp = IC_cpr1D_oneslice./RMStmp;
        CSP_tmp = (CSP_tmp./CSP_tmp(:,:,1)).*abs(CSP_tmp(:,:,1));
        
        CSP_tmp_pc(:,:,k1,k2,:) = CSP_tmp;
        IC_tmp_pc(:,:,k1,k2,:) = IC_cpr1D_oneslice;
    end
end

RMS_tmp_pc = zeros(Ny_MB,Nx_MB,length(pc1),size(pc2,1));
for k1 = 1 : length(pc1)
    for k2 = 1 : size(pc2,1)
        RMS_tmp_pc(:,:,k1,k2) = RootMeanSquare(squeeze(IC_tmp_pc(:,:,k1,k2,:)),1,Ncoil_SB,'Image');
    end
end

Mask_CSP_cpr1D = Mask_csp;
Mask_bg_cpr1D = Mask_bg;

Difftmp_csp = zeros(length(pc1),size(pc2,1),Ncoil_MB);
Difftmp_bg = zeros(length(pc1),size(pc2,1),Ncoil_MB);
for c = 1 : Ncoil_MB
    q1_s1 = CSP_tmp_pc(:,:,:,:,c).*Mask_CSP_cpr1D;
    q2_s1 = squeeze(CSP_MB(:,:,choosethisslice,c)).*Mask_CSP_cpr1D;
    
    Difftmp_bg = RMS_tmp_pc.*Mask_bg_cpr1D;

    Diff_s1 = abs(abs(q1_s1)-abs(repmat(q2_s1,1,1,length(pc1),size(pc2,1))));
    Difftmp_csp(:,:,c) = squeeze(sum(sum(Diff_s1,1),2));
end
tmp = squeeze(sqrt(mean(abs(Difftmp_csp).^2,3))); % Combine coils;
[p1_csp_zz,p2_csp_zz] = find(tmp == min(tmp(:)));

tmp1 = squeeze(sum(sum(Difftmp_bg,1),2));
[p1_bg_zz,p2_bg_zz] = find(tmp1 == min(tmp1(:)));

pc1_zz_cpr1D = pc1(p1_csp_zz).*weight_csp + pc1(p1_bg_zz).*(1 - weight_csp);
pc2_zz_cpr1D = pc2(p2_csp_zz,:).*weight_csp + pc2(p2_bg_zz,:).*(1 - weight_csp);
IC_cpr1D_oneslice = zeros(Ny_MB,Nx_MB,Ncoil_MB);
kspace_cpr1D_oneslice = zeros(Ny_MB,Nx_MB,Ncoil_MB);

for c = 1 : Ncoil_MB    
    for y = 1 : nshot
        arrtmp = kspace_MB_2shot(y+nshot:ns4:end,:,choosethisslice,c)*pc1_zz_cpr1D;
        
        ktmp = fftshift(fft(fftshift(arrtmp,2),[],2),2);
        pc2tmp = pc2_zz_cpr1D;
        pc2tmp = repmat(pc2tmp,[Ny_MB/ns4 1]);
        kcorr = ktmp.*pc2tmp;
        kspace_cpr1D_oneslice(y+nshot:ns4:end,:,c) = ifftshift(ifft(ifftshift(kcorr,2),[],2),2);
        kspace_cpr1D_oneslice(y:ns4:end,:,c) = kspace_MB_2shot(y:ns4:end,:,choosethisslice,c);
    end
    IC_cpr1D_oneslice(:,:,c) = fftshift(fft2(fftshift(kspace_cpr1D_oneslice(:,:,c))));
end

RMS_cpr1D_oneslice = RootMeanSquare(squeeze(IC_cpr1D_oneslice),1,Ncoil_SB,'Image');
clear IC_tmp_pc;

%% Run CPR on one slice: cpr1DB non-linear
SaveFactor3 = 0.1;
stepN3 = 31;
coeff3step = SaveFactor3*2*pi/stepN3;
coeff3list = -(SaveFactor3*pi) +  coeff3step: coeff3step : (SaveFactor3*pi);
pc3 = exp(1i*coeff3list);

kspace_tmp = zeros(Ny_MB,Nx_MB,length(pc3),Ncoil_MB);
IC_tmp_pc3 = zeros(Ny_MB,Nx_MB,length(pc3),Ncoil_MB);
CSP_tmp_pc3 = zeros(Ny_MB,Nx_MB,length(pc3),Ncoil_MB);
RMS_tmp_pc3 = zeros(Ny_MB,Nx_MB,length(pc3));

for k1 = 1 : length(pc3)
    for c = 1 : Ncoil_MB
        kspace_tmp(:,:,k1,c) = kspace_cpr1D_oneslice(:,:,c);
        for y = 1 : nshot
            kspace_tmp(y+nshot:ns4:end,:,k1,c) = kspace_cpr1D_oneslice(y+nshot:ns4:end,:,c).*pc3(k1);
            kspace_tmp(y:ns4:end,:,k1,c) = kspace_cpr1D_oneslice(y:ns4:end,:,c);
        end
        
        IC_tmp_pc3(:,:,k1,c) = fftshift(fft2(fftshift(kspace_tmp(:,:,k1,c))));
    end
    RMS_tmp_pc3(:,:,k1) = RootMeanSquare(squeeze(IC_tmp_pc3(:,:,k1,:)),1,Ncoil_MB,'Image');
    CSP_tmp_pc3(:,:,k1,:) = squeeze(IC_tmp_pc3(:,:,k1,:))./RMS_tmp_pc3(:,:,k1);
    CSP_tmp_pc3(:,:,k1,:) = (CSP_tmp_pc3(:,:,k1,:)./CSP_tmp_pc3(:,:,k1,1)).*abs(CSP_tmp_pc3(:,:,k1,1));
end

pc3_csp_zz = zeros(1,Nx_MB);
pc3_bg_zz = zeros(1,Nx_MB);
p3_csp_zz = zeros(1,Nx_MB);
p3_bg_zz = zeros(1,Nx_MB);

Diffsum_pc3 = zeros(Nx_MB,length(pc3));
Diffsum_bg_pc3 = zeros(Nx_MB,length(pc3));

Mask_CSP_cpr1DB = Mask_csp;
Mask_bg_cpr1DB = Mask_bg;

for x = xstart_bg : xend_bg
    Difftmp_csp = zeros(length(pc3),Ncoil_MB);
    if sum(Mask_CSP_cpr1DB(:,x)) > 0
        for c = 1 : Ncoil_MB
            q1_s1 = squeeze(CSP_tmp_pc3(:,x,:,c)).*Mask_CSP_cpr1DB(:,x);
            q2_s1 = squeeze(CSP_MB(:,x,choosethisslice,c)).*Mask_CSP_cpr1DB(:,x);
            
            Diff_s1 = abs(abs(q1_s1)-abs(repmat(q2_s1,1,length(pc3))));
            Difftmp_csp(:,c) = squeeze(sum(Diff_s1,1));
        end
        
        tmp = squeeze(sqrt(mean(abs(Difftmp_csp).^2,2))); % Combine coils;
        p3_csp_zz(1,x) = find(tmp == min(tmp(:)));
        Diffsum_pc3(x,:) = tmp;
        pc3_csp_zz(1,x) = pc3(p3_csp_zz(x));
        
        pc3_bg_zz(1,x) = pc3(p3_csp_zz(x));        
    else        
        tmp1 = sum(squeeze(RMS_tmp_pc3(:,x,:)).*Mask_bg_cpr1DB(:,x),1);
        Diffsum_bg_pc3(x,:) = tmp1;
        p3_bg_zz(1,x) = find(tmp1 == min(tmp1(:)));
        pc3_bg_zz(1,x) = pc3(p3_bg_zz(x));
        pc3_csp_zz(1,x) = pc3(p3_bg_zz(x));        
    end     
end

pc3_zz = pc3_csp_zz.*weight_csp + pc3_bg_zz.*(1 - weight_csp);
outliner = isoutlier(angle(pc3_zz(xstart_bg : xend_bg)),'quartiles');
index = xstart_bg : xend_bg;
index = index(outliner~=1);

% Unwrap and fitting
savemepcangle = unwrap(angle(pc3_zz(outliner~=1)));
xlist = index;

p = polyfit(xlist,savemepcangle,2);
pc3_zz_cpr1DB_pro = exp(1i*polyval(p,1:Nx_MB));

% Correcttion after fitting
kspace_cpr1DB_oneslice = zeros(size(kspace_cpr1D_oneslice));
IC_cpr1DB_oneslice = zeros(size(kspace_cpr1DB_oneslice));
for c = 1 : Ncoil_MB
    kspace_cpr1DB_oneslice(:,:,c) = kspace_cpr1D_oneslice(:,:,c);
    for y = 1 : nshot
        arrtmp = kspace_cpr1D_oneslice(y+nshot:ns4:end,:,c);
        
        ktmp = fftshift(fft(fftshift(arrtmp,2),[],2),2);
        kcorr = ktmp.*pc3_zz_cpr1DB_pro;
        kspace_cpr1DB_oneslice(y+nshot:ns4:end,:,c) = ifftshift(ifft(ifftshift(kcorr,2),[],2),2);
        kspace_cpr1DB_oneslice(y:ns4:end,:,c) = kspace_cpr1D_oneslice(y:ns4:end,:,c);
    end
    
    IC_cpr1DB_oneslice(:,:,c) = fftshift(fft2(fftshift(kspace_cpr1DB_oneslice(:,:,c))));
end
RMS_cpr1DB_oneslice = RootMeanSquare(IC_cpr1DB_oneslice,1,Ncoil_MB,'Image');
clear IC_tmp_pc3;

 %% 2D correction
PhaseMap_1D = repmat(angle(pc1_zz_cpr1D.*pc2_zz_cpr1D.*pc3_zz_cpr1DB_pro),Ny_MB,1);
[Phase_coeff,PhaseMap_pc] = cpr2D_coils_mb2_0930_2s(squeeze(kspace_MB_2shot(:,:,choosethisslice,:)),squeeze(Image_SB_ref(:,:,MB_SliceNum(:,choosethisslice),:)),squeeze(CSP_MB(:,:,choosethisslice,:)),squeeze(CSP_SB(:,:,MB_SliceNum(:,choosethisslice),:)),PhaseMap_1D,FOVshift,nshot,weight_csp);
PhaseMap_1D_MB = PhaseMap_1D;
PhaseMap_2D_MB = PhaseMap_pc;
save([FolderPath '/PhaseMap_1D_MB.mat'],'PhaseMap_1D_MB');
save([FolderPath '/PhaseMap_2D_MB.mat'],'PhaseMap_2D_MB');

%% Slice-dependent Correction (will take longer time)
% [CorrImg_SD,PhaseMap_SD,Phase_coeff_SD] = cpr2D_coils_mb2_SD_0930_2s(squeeze(kspace_MB_2shot(:,:,choosethisslice,:)),squeeze(Image_SB_ref(:,:,MB_SliceNum(:,choosethisslice),:)),squeeze(CSP_SB(:,:,MB_SliceNum(:,choosethisslice),:)),Phase_coeff,PhaseMap_1D,FOVshift,nshot);
% PhaseMap_SD_MB = PhaseMap_SD;
% save([FolderPath '/PhaseMap_SD_MB.mat'],'PhaseMap_SD_MB');
% PhaseMap = PhaseMap_SD_MB;

clear kspace_MB_2shot;

PhaseMap = repmat(PhaseMap_2D_MB,[1 1 2]);
Corr_MB_2D_MUSE_2s;

