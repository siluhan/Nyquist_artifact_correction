function [Mask_all,Mask_threshold_SB] = Mask_Generate(RMS_all,StartThreshold,EndThreshold)


[Ny_SB,Nx_SB,Nz_SB] = size(RMS_all);
for z = 1 : Nz_SB
    hanning20 = hanning(floor(Nx_SB/8))*hanning(floor(Nx_SB/8))';
    tmp1 = RMS_all(:,:,z);  % Select one image
    tmp2 = conv2(tmp1,hanning20,'same'); % Image filtered with a hanning filter
    Mask_s0 = tmp2 > StartThreshold*max(tmp2(:));
    
    Mask_s0 = imfill(Mask_s0,'holes');
    A = regionprops(Mask_s0);
    Area = zeros(length(A),1);
    for k = 1 : length(A)
        Area(k) = A(k).Area;
    end
    Area = sort(Area,'descend');
    if length(A) > 1
        Mask_s0 = bwareaopen(Mask_s0,Area(2)+1);
    end
    
    Mask_threshold_SB = StartThreshold;
    while sum(sum(~Mask_s0(1:Ny_SB/2,:),1) < ceil(Ny_SB*0.1/2)) > 0 || sum(sum(~Mask_s0(Ny_SB/2+1:end,:),1) < ceil(Ny_SB*0.1/2)) > 0 || sum(sum(~Mask_s0,1) < ceil(Ny_SB*0.1)) > 0 %To make every column has at least 5 pixels for bg at both top and bottom or total 10 pixels
        Mask_threshold_SB = Mask_threshold_SB + 0.05;
        Mask_s0 = tmp2 > Mask_threshold_SB*max(tmp2(:));
        
        Mask_s0 = imfill(Mask_s0,'holes');
        A = regionprops(Mask_s0);
        Area = zeros(length(A),1);
        for k = 1 : length(A)
            Area(k) = A(k).Area;
        end
        Area = sort(Area,'descend');
        if length(A) > 1
            Mask_s0 = bwareaopen(Mask_s0,Area(2)+1);
        end
        if Mask_threshold_SB > EndThreshold
            break;
        end
    end
    
    Mask_threshold_SB_MUSE = Mask_threshold_SB + 0.05;
    
    Mask_s0 = tmp2 > Mask_threshold_SB_MUSE*max(tmp2(:));
    
    Mask_s0 = imfill(Mask_s0,'holes');
    A = regionprops(Mask_s0);
    Area = zeros(length(A),1);
    for k = 1 : length(A)
        Area(k) = A(k).Area;
    end
    Area = sort(Area,'descend');
    if length(A) > 1
        Mask_s0 = bwareaopen(Mask_s0,Area(2)+1);
    end
    Mask_all(:,:,z) = Mask_s0;
end
