function denseSIFT = denseSIFT(file, gridSpacing)
    fprintf('Building Sift Descriptors\n\n');

    I=imread(file); % read the image
    if ndims(I) == 3
        I = im2double(rgb2gray(I));
    else
        I = im2double(I);
    end

    [hgt, wid] = size(I);


    %% make grid (coordinates of upper left patch corners)
    %remX = mod(wid-patchSize,gridSpacing);% the right edge
    %offsetX = floor(remX/2)+1;
    %remY = mod(hgt-patchSize,gridSpacing);
    %offsetY = floor(remY/2)+1;
    %[gridX,gridY] = meshgrid(offsetX:gridSpacing:wid-patchSize+1,offsetY:gridSpacing:hgt-patchSize+1);

    %% make grid edited
    remX = mod(wid,gridSpacing);% the right edge
    offsetX = floor(remX/2)+1;
    remY = mod(hgt,gridSpacing);
    offsetY = floor(remY/2)+1;
    [gridX,gridY] = meshgrid(offsetX:gridSpacing:wid,offsetY:gridSpacing:hgt);

    bin = [];
    for i = 1:size(gridX,1)
        for j = 1:size(gridX,2)
            bin = [bin;[gridX(i,j),gridY(i,j)]];
        end
    end

    orient = buildDescriptor(I, bin, gridSpacing);
    %writeLineDesc(orient)
    showFigure(I, gridSpacing,orient)
    %gridImage(I, gridSpacing)
    
    denseSIFT = lineDesc(orient);

function buildDescriptor = buildDescriptor(I, bin, gridSpacing)
    
    %for each of the keypoints, calculate the orientation, then rotate
    %the keypoint descriptor accordingly. Finally, calculate the keypoint
    %descriptor. 
    cont = 0;
    qtyKeypoints = size(bin,1);
    
    kptDescriptors = repmat(struct('kptX',0,'kptY',0,'kptDescriptor',zeros(4,4,8)),qtyKeypoints,1);

    %descriptor window size 
    DESC_WIN_SIZE = 16;

    %filter for calculating diffX: 
    filterDiffX = [0 0 0; -1 0 1; 0 0 0];

    diffXMat = imfilter(I, filterDiffX);

    %filter for calculating diffY: 
    filterDiffY = [0 1 0; 0 0 0; 0 -1 0];

    diffYMat = imfilter(I, filterDiffY);

    %get the magnitude operating directly on matrixes: 
    magnMat = sqrt(diffXMat.*diffXMat + diffYMat.*diffYMat);

    %do similar thing for orientation
    orientMat = atan2(diffYMat, diffXMat);

    for binCenter = 1:size(bin,1)
        
        cont = cont+1;
        %the gaussian weights for the window 
        gaussWeight = getGaussWeights(DESC_WIN_SIZE, DESC_WIN_SIZE/2);

        %%%%%%%%%gets the coordinates image (of magnitudes)%%%%%%
        r = gridSpacing/2;
        row = bin(binCenter,1) + r - 0.5;
        col = bin(binCenter,2) + r - 0.5;
        
        kptDescriptor = zeros(128,1);
        

        %the window is 16 x 16 pixels in the keypoint level
        for x = 0:DESC_WIN_SIZE-1

            for y = 0:DESC_WIN_SIZE-1
                %first identify subregion I am in 
                subregAxisX = floor(x/4);
                subregAxisY = floor(y/4);

                yCoord = row + y - DESC_WIN_SIZE/2; 
                xCoord = col + x - DESC_WIN_SIZE/2; 
                yCoord = round(yCoord); 
                xCoord = round(xCoord); 
                %get the magnitude

                if(yCoord>0&&xCoord>0&&yCoord<=size(magnMat,1) && xCoord<=size(magnMat,2)) 

                    magn = magnMat(yCoord,xCoord); 

                    %multiply the magnitude by gaussian weight 
                    magn = magn*gaussWeight(y+1,x+1); 

                    orientation = orientMat(yCoord,xCoord);
                    orientation = orientation + pi;
                    %calculate the respective bucket

                    bucket = (orientation)*(180/pi);
                    bucket = ceil(bucket/45);

                    kptDescriptor((subregAxisY+subregAxisX*4)*8 + bucket) = ...
                        kptDescriptor((subregAxisY+subregAxisX*4)*8 + bucket) + magn;
                end
            end
        end
        
        %normalize the vector
        sqKptDescriptor = kptDescriptor.^2; 
        sumSqKptDescriptor = sum(sqKptDescriptor);
        dem = sqrt(sumSqKptDescriptor); 
        kptDescriptor = kptDescriptor./dem;
        
        %fprintf('%i>',binCenter);
        for i=1:size(kptDescriptor,1)
            %fprintf(' %i',kptDescriptor(i));
            if isnan(kptDescriptor(i))
                kptDescriptor(i) = 0;
            end
        end
        %fprintf('\n');
        
        %threshold 
        kptDescriptor(find(kptDescriptor>0.2))=0.2;
        
        %Renormalizing again, as stated in 6.1 of [1]
        sqKptDescriptor = kptDescriptor.^2; 
        sumSqKptDescriptor = sum(sqKptDescriptor);
        dem = sqrt(sumSqKptDescriptor); 
        kptDescriptor = kptDescriptor./dem;
        
        for i=1:size(kptDescriptor,1)
            %fprintf(' %i',kptDescriptor(i));
            if isnan(kptDescriptor(i))
                kptDescriptor(i) = 0;
            end
        end
        
        %keypointDescriptors{octave}{kptLayer}{rowKpt(keypoint)}{rowKpt(keypoint)} = kptDescriptor;
        kptDescriptors(cont) = struct('kptX',row,'kptY',col, ...
            'kptDescriptor',kptDescriptor);

    end
    %save('defineOrient.mat','magnMat','orientMat');

    buildDescriptor = kptDescriptors;
    
function getGaussWeights = getGaussWeights(windowSize, sigma)
    k = fspecial('Gaussian', [windowSize windowSize], sigma);
    k = k.*(1/max(max(k))); 
    getGaussWeights = k;

function showFigure(I, gridSpacing,desc)
    [rows, columns, ~] = size(I);
    I = 255 * ones(rows, columns, 'uint8');
    iptsetpref('ImshowBorder','tight');
    imshow(I);
    axis on;
    hold on; 
    for row = 1 : gridSpacing : rows
        line([1, columns], [row, row], 'Color', 'r', 'LineWidth', 0.5);
    end
    for col = 1 : gridSpacing : columns
        line([col, col], [1, rows], 'Color', 'r', 'LineWidth', 0.5);
    end
    
    cont=1;
    for i=1:length(desc)
        kpt=extractfield(desc(i),'kptDescriptor');
        x=extractfield(desc(i),'kptX');
        y=extractfield(desc(i),'kptY');
        %fprintf('%i = ',i)
        temp = createLocs(x,y,kpt);
        for j=1:length(temp)
            %grayImage(cent(1)+(-1:1), cent(2)+(-1:1)) = grayImage(cent(1)+(-1:1), cent(2)+(-1:1));
            %temp(j,1:end)
            locs(cont,1:size(temp,2))=temp(j,1:end);
            cont=cont+1;
        end
        storeDesc(i,:) = kpt(:);
    end
    
    for i=1:length(locs)
        % Draw an arrow, each line transformed according to keypoint parameters.
        TransformLine(locs(i,:), 0.0, 0.0, 1.0, 0.0);
        TransformLine(locs(i,:), 0.85, 0.1, 1.0, 0.0);
        TransformLine(locs(i,:), 0.85, -0.1, 1.0, 0.0);
    end
    hold off;
    saveas(gcf,'imagegridwithoutimage.png')

function writeLineDesc(desc)
    %
    for i=1:length(desc)
        fprintf('%i ',i);
        kpt=extractfield(desc(i),'kptDescriptor');
        for j=1:length(kpt)
            fprintf('%i:%5.5f ',j,kpt(j));
        end
        fprintf('\n');
    end

function res = lineDesc(desc)
    %
    for i=1:length(desc)
        kpt=extractfield(desc(i),'kptDescriptor');
        storeDesc(i,:) = kpt(:);
    end
    
    res = reshape(storeDesc.',1,[]);
    res = res.';

function gridImage(image, gridSpacing)
    %
    folder = 'C:\Users\Alex Winter\Pictures\SplitImage\';
    [m,n] = size(image);
    blocks = cell(m/gridSpacing,n/gridSpacing);
    count = 0;
    counti = 0;
    for i = 1:16:m-15
       counti = counti + 1;
       countj = 0;
       for j = 1:16:n-15
            countj = countj + 1;
            blocks{counti,countj} = image(i:i+15,j:j+15);
            count = count + 1;
            %fprintf('%i Loop i=%i,j=%i\n',count,i,j);
            
            %export images
            filename = sprintf('Split Image %i.jpg', count);
            imwrite(blocks{counti,countj}, fullfile(folder, filename));
       end
    end

function createLocs = createLocs(x,y,desc)
    %
    ori=-(0:7)*pi/4;
    for i=1:length(desc)
        row = ceil(ceil(i/8)/4);
        col = mod(ceil(i/8)-1,4)+1;
        %fprintf('%i, ',ceil(i/8))
        cont = i-(8*(ceil(i/8)-1));
        %fprintf('%i %5.2f %5.2f %i\n',i,x+(col*8)-20,y+(row*8)-20,cont);
        loc(i,:)=[x+(col*8)-20,y+(row*8)-20,desc(i),ori(cont)];
    end
    %fprintf('\n');
    createLocs = loc;
    
function TransformLine(keypoint, x1, y1, x2, y2)
    
    
    % The scaling of the unit length arrow is set to approximately the radius
    %   of the region used to compute the keypoint descriptor.
    len = 6 * keypoint(3);

    % Rotate the keypoints by 'ori' = keypoint(4)
    s = sin(keypoint(4));
    c = cos(keypoint(4));

    % Apply transform
    r1 = keypoint(1) - len * (c * y1 + s * x1);
    c1 = keypoint(2) + len * (- s * y1 + c * x1);
    r2 = keypoint(1) - len * (c * y2 + s * x2);
    c2 = keypoint(2) + len * (- s * y2 + c * x2);

    grid on
    grid minor

    line([c1 c2], [r1 r2], 'Color', 'c');