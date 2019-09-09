%function to make gaussian blur of the image, the kernel size is calculated automatically 
function gaussianBlur=gaussianBlur(image, sigma)


    kernelSize = round(sigma*3 - 1); 
    kernelSize
    if(kernelSize<1)
        kernelSize = 1; 
    end 
	kernel = fspecial('gaussian', [kernelSize kernelSize], sigma);
    
    kernel
    
	convImage = imfilter(image,kernel,'replicate');
%	convImage = imfilter(image,kernel);

	gaussianBlur = convImage; 
end