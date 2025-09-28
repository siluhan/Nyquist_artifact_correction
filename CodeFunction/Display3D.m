%% Display 3D images in 2D figure
function Display3D(Image,ColumnNum,Order)
% Image is a 3D matrix. ColumeNum is defined as the number of images
% in colume to display. Order is the increment of image number to 
% display. The number of images in row is automatically determined.
% Contrast is a matrix containing [cmin cmax].

% Edited on Jan 12, 2018

ImageSize = size(Image);
if length(ImageSize) == 3
    RowNum = ceil(ImageSize(3)/Order/ColumnNum);
    DisplayImage = zeros(RowNum*ImageSize(1),ColumnNum*ImageSize(2));
    n = 0;
    for k = 1 : Order : ImageSize(3)
        n = n + 1;
        Row = ceil(n/ColumnNum);
        Colume = n - (Row - 1)*ColumnNum;
        DisplayImage((Row - 1)*ImageSize(1) + 1 : Row*ImageSize(1), (Colume - 1)...
            *ImageSize(2) + 1 : Colume*ImageSize(2)) = Image(:,:,k);
    end
    imshow(DisplayImage,[]);
elseif length(ImageSize) == 2
    imshow(Image,[])
end



    
    
    