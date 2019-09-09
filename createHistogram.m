function createHistogram(file)
    %
    load('defineOrient.mat')
    I=imread(file); % read the image
    if ndims(I) == 3
        I = im2double(rgb2gray(I));
    else
        I = im2double(I);
    end
    
    
    imshow(I);
    axis on;
    hold on;
    for i=1:size(I,1)
        for j=1:size(I,2)
            locs = [i,j,magnMat(i,j),orientMat(i,j)];
            % Draw an arrow, each line transformed according to keypoint parameters.
            TransformLine(locs, 0.0, 0.0, 1.0, 0.0);
            TransformLine(locs, 0.85, 0.1, 1.0, 0.0);
            TransformLine(locs, 0.85, -0.1, 1.0, 0.0);
        end
    end
    hold off;

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