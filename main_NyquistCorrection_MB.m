clear variables;

currentFolder = pwd;
addpath([currentFolder '/CodeFunction']);

%% Load GRE kspace data
FolderPath = [currentFolder '/example_data'];  
load([FolderPath '/kspace_gre_demo_MB.mat' ]);

kspace_gre = kspace;
[Ny_SB,Nx_SB,Nz_SB,Ncoil_SB,Nt_SB] = size(kspace_gre);
clear kspace;

Image_gre = zeros(size(kspace_gre));

for z = 1 : Nz_SB
    for c = 1 : Ncoil_SB
        Image_gre(:,:,z,c) = fftshift(fft2(fftshift(kspace_gre(:,:,z,c))));
    end
    
end

RMS_gre = RootMeanSquare(Image_gre,Nz_SB,Ncoil_SB,'Image');
figure;Display3D(imrotate(squeeze(RMS_gre),0),13,1);title('GRE')


%% Read EPI data: 2-shot EPI
FolderPath = [currentFolder '/example_data'];  
load([FolderPath '/kspace_epi_2shot_MB.mat' ]);

% kspace_epi_mb = kspace(:,:,:,:,6:end); % skip first 5 tpns
kspace_epi_mb = kspace;
clear kspace;

[~,~,Nz_MB,Ncoil_MB,Nt_MB] = size(kspace_epi_mb);
Image_epi = zeros(size(kspace_epi_mb));
RMS_epi = zeros(size(kspace_epi_mb,[1 2 3 5]));
for t = 1 : Nt_MB
    for z = 1 : Nz_MB
        for c = 1 : Ncoil_MB
            Image_epi(:,:,z,c,t) = fftshift(fft2(fftshift(kspace_epi_mb(:,:,z,c,t))));
        end
    end

    RMS_epi(:,:,:,t) = RootMeanSquare(squeeze(Image_epi(:,:,:,:,t)),Nz_MB,Ncoil_MB,'Image'); 

end
figure;Display3D(imrotate(RMS_epi(:,:,:,Nt_MB),0),1,1); title(['2-shot EPI, tpn = ' num2str(Nt_MB) ])

%% coil sensitivity profile
CSP_SB = Image_gre./RMS_gre;
CSP_SB = (CSP_SB./CSP_SB(:,:,:,1)).*abs(CSP_SB(:,:,:,1));

zz_mb = [1:Nz_SB/2;Nz_SB/2+1:Nz_SB];
Phase_mb2_s1 = repmat(repmat([0;pi],1,Nx_SB),Ny_SB/2,1); % 1/2 FOV shift
kspace_gre_mb = zeros(Ny_SB,Nx_SB,Nz_MB,Ncoil_SB);
Image_gre_mb = zeros(Ny_SB,Nx_SB,Nz_MB,Ncoil_SB);
for k1 = 1 : Nz_MB
    for k2 = 1 : Ncoil_SB
        kspace_s1 = kspace_gre(:,:,zz_mb(2,k1),k2).*exp(1i*Phase_mb2_s1);
        
        kspace_gre_mb(:,:,k1,k2) = kspace_gre(:,:,zz_mb(1,k1),k2) + kspace_s1;
        
        Image_gre_mb(:,:,k1,k2) = fftshift(fft2(fftshift(kspace_gre_mb(:,:,k1,k2))));
    end
end
RMS_gre_mb = RootMeanSquare(Image_gre_mb,Nz_MB,Ncoil_SB,'Image');

CSP_MB = Image_gre_mb./RMS_gre_mb;
CSP_MB = (CSP_MB./CSP_MB(:,:,:,1)).*abs(CSP_MB(:,:,:,1));

%% 2-shot MB-EPI correction
Image_SB_ref = Image_gre;
RMS_SB_ref = RMS_gre;
kspace_MB_tpn = kspace_epi_mb;

[Ny_MB,Nx_MB,Nz_MB,Ncoil_MB,Nt_MB] = size(kspace_MB_tpn);
kdataoriginal = mean(kspace_MB_tpn,5);
rms = RootMeanSquare(kdataoriginal,Nz_MB,Ncoil_SB,'kspace');

nshot = 2;
ns4 = nshot*2;
FOVshift = 1/2; 
weight_csp = 1;   % determine phase value purely based on coil sensitivity profile

Corr_MB_2shot_csp_MUSE; % Perform 2D phase correction

% Display results
figure;subplot(1,2,1);Display3D(abs(RMS_epi(:,:,:,5)),1,1);title('original 2-shot EPI')
subplot(1,2,2);Display3D(abs(CorrImg_MB_2shot_2D(:,:,:,5)),2,1);title('corrected 2-shot EPI')
figure;subplot(1,2,1);Display3D(abs(RMS_epi(:,:,:,5)),1,1);title('original 2-shot EPI');clim([0 max(RMS_epi(:))/6]);
subplot(1,2,2);Display3D(abs(CorrImg_MB_2shot_2D(:,:,:,5)),2,1);title('corrected 2-shot EPI');clim([0 max(RMS_epi(:))/6]);

