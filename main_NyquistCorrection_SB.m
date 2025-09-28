clear variables;

currentFolder = pwd;
addpath([currentFolder '/CodeFunction']);

%% Load GRE kspace data
FolderPath = [currentFolder '/example_data'];  
load([FolderPath '/kspace_gre_demo_SB.mat' ]);

[~,~,Nz_SB,Ncoil_SB] = size(kspace);
kspace_gre = kspace;
Image_gre = zeros(size(kspace_gre));

for z = 1 : Nz_SB
    for c = 1 : Ncoil_SB
        Image_gre(:,:,z,c) = fftshift(fft2(fftshift(kspace_gre(:,:,z,c))));
    end    
end

RMS_gre = RootMeanSquare(Image_gre,Nz_SB,Ncoil_SB,'Image');
figure;Display3D(imrotate(squeeze(RMS_gre),0),2,1);title('GRE')

CSP_SB = Image_gre./RMS_gre;
CSP_SB = (CSP_SB./CSP_SB(:,:,:,1)).*abs(CSP_SB(:,:,:,1));

%% Read EPI data: 2-shot EPI
FolderPath = [currentFolder '/example_data'];  
load([FolderPath '/kspace_epi_2shot_SB.mat' ]);

% kspace_epi = kspace(:,:,:,:,6:end); % skip first 5 tpns
kspace_epi = kspace;
[Ny_SB,Nx_SB,Nz_SB,Ncoil_SB,Nt_SB] = size(kspace_epi);
clear kspace;

Image_epi = zeros(size(kspace_epi));
RMS_epi = zeros(size(kspace_epi,[1 2 3 5]));
for t = 1 : Nt_SB
    for z = 1 : Nz_SB
        for c = 1 : Ncoil_SB
            Image_epi(:,:,z,c,t) = fftshift(fft2(fftshift(kspace_epi(:,:,z,c,t))));
        end
    end
    RMS_epi(:,:,:,t) = RootMeanSquare(squeeze(Image_epi(:,:,:,:,t)),Nz_SB,Ncoil_SB,'Image'); 
end
figure;Display3D(imrotate(RMS_epi(:,:,:,Nt_SB),0),2,1);
title(['2-shot EPI, tpn = ' num2str(Nt_SB) ]);

%% 2-shot EPI correction
Image_SB_ref = Image_gre;
RMS_SB_ref = RMS_gre;
kspace_SB_tpn = kspace_epi;

kdataoriginal = mean(kspace_SB_tpn,5);

nshot = 2;
ns4 = nshot*2;
weight_csp = 1; % determine phase value purely based on coil sensitivity profile

Corr_SB_2shot_csp_MUSE; % Perform 2D phase correction

% Display results
choosethisslice = ceil(Nz_SB/2);
figure;subplot(1,2,1);imshow(abs(RMS_epi(:,:,choosethisslice,5)),[]);title('original 2-shot EPI')
subplot(1,2,2);imshow(abs(CorrImg_SB_2shot_2D(:,:,choosethisslice,5)),[]);title('corrected 2-shot EPI')
figure;subplot(1,2,1);imshow(abs(RMS_epi(:,:,choosethisslice,5)),[]);title('original 2-shot EPI');clim([0 max(RMS_epi(:))/6]);
subplot(1,2,2);imshow(abs(CorrImg_SB_2shot_2D(:,:,choosethisslice,5)),[]);title('corrected 2-shot EPI');clim([0 max(RMS_epi(:))/6]);

