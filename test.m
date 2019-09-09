%image = imread('D:\OneDrive\Document\SIFT(Scale-Invariant Feature Transform)\Code\denseSIFT\img\A014.jpg');
%image = rgb2gray(image);
%blurImage = gaussianBlur(image, 1.8);
%figure
%subplot(1,2,1);
%imshow(image);
%subplot(1,2,2);
%imshow(blurImage);
fprintf('%i %i %i\n')
for i=1:128
    row = ceil(ceil(i/8)/4);
    col = mod(ceil(i/8)-1,4)+1;
    fprintf('%i %i %i %i\n',i,ceil(i/8),(col*8)-20,(row*8)-20);
end