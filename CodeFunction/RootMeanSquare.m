%% Generate RMS Image
function RMSsum = RootMeanSquare(Input,Nz,Ncoil,DataType)
% The input can be image or kspace. 
% The size of input image/kspace is [Nx,Ny,Nz,Ncoil]
% DataType can be: Image or kspace


Nx = size(Input,1);
Ny = size(Input,2);

if Nz > 1
    ImageSum = zeros(Nx,Ny,Nz,Ncoil);
    RMSsum = zeros(Nx,Ny,Nz);
    
    if strcmp(DataType,'Image')
        for k1 = 1 : Nz
            RMS = sqrt(mean(abs(squeeze(Input(:,:,k1,:))).^2,3));
            RMSsum(:,:,k1) = RMS;
        end
        
    elseif strcmp(DataType,'kspace')
        for k1 = 1 : Nz
            for k2 = 1 : Ncoil
                Image = fftshift(fft2(fftshift(Input(:,:,k1,k2))));
                ImageSum(:,:,k1,k2) = Image;
            end
            RMS = sqrt(mean(abs(squeeze(ImageSum(:,:,k1,:))).^2,3));
            RMSsum(:,:,k1) = RMS;
        end
    end
    
elseif Nz == 1
       ImageSum = zeros(Nx,Ny,Ncoil);
       RMSsum = zeros(Nx,Ny);
    
    if strcmp(DataType,'Image')

            RMS = sqrt(mean(abs(squeeze(Input)).^2,3));
            RMSsum = RMS;

        
    elseif strcmp(DataType,'kspace')
            for k2 = 1 : Ncoil
                Image = fftshift(fft2(fftshift(Input(:,:,k2))));
                ImageSum(:,:,k2) = Image;
            end
            RMS = sqrt(mean(abs(squeeze(ImageSum)).^2,3));
            RMSsum = RMS;
    end
    
end





