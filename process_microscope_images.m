%
% Syntax:  [output] = function_name(inputs)
%
% Inputs:
%   1.
%
% Outputs:
%	1.
%
% Example:
%
%
% Other m-files required:
%
% Subfunctions:
%
% MAT-files required:
%
% See also:
%
% Author: Brandon patterson
% email address: i.am.brandon.patterson@gmail.com
%
% Copyright (C) 2018, Brandon Patterson
%
% Last revision: 02-07-2018

%% Sort image by color


%% Sort image by color
if true
    clc; clear all; close all;
    
    % Read the image file
    ImDir='e:/brandon/microscope_images/';
    ImFile = 'LTA1_x40_edge_of_hem.jpg';
    Im = imread(sprintf('%s%s',ImDir,ImFile));
    
    % Convert RGB to L*a*b
    lab_im = rgb2lab(Im);
    
    %Calculate the mean 'a*' and 'b*' value for each area that you extracted with roipoly. These values serve as your color markers in 'a*b*' space.
    a = lab_im(:,:,2);
    b = lab_im(:,:,3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET COLORS %%%%%%%%%%%%%%%%%%%%
    
    % Define a function to grab a square for each color
%     mysquare = @(x,y,s) ([x+floor(s/2).*[-1;-1;+1;+1;-1], y+floor(s/2).*[-1;+1;+1;-1;-1]]);
%     
%     % Rgioncoordinates for areas of each color (this is done manually right now)
%     white_region = mysquare(762,503,14);
%     dpurple_region = mysquare(849,565,10);
%     ddpurple_region = mysquare(875,947,10);
%     lpurple_region = mysquare(790,602,8);
%     transparent_region = mysquare(1918,709,10);
%     red_region = mysquare(572,836,10);
%     
%     region_coordinates(:,:,1) = white_region;
%     region_coordinates(:,:,2) = dpurple_region;
%     region_coordinates(:,:,3) = lpurple_region;
%     region_coordinates(:,:,4) = transparent_region;
    % region_coordinates(:,:,5) = ddpurple_region; %UNCOMMENTING THIS TO ACCOUNT FOR DARK PURPLE BREAKS EVERYTHING
    
%     ncolors = 4;
%     sample_regions = false([size(Im,1) size(Im,2) ncolors]);
    
%     for count = 1:ncolors
%         sample_regions(:,:,count) = roipoly(Im,region_coordinates(:,1,count), region_coordinates(:,2,count));
%     end
%     
    
    
    if false
        color_markers = zeros([ncolors, 2]);
        for count = 1:ncolors
            color_markers(count,1) = mean2(a(sample_regions(:,:,count)));
            color_markers(count,2) = mean2(b(sample_regions(:,:,count)));
        end
        
        %For example, the average color of the dark purple sample region in 'a*b*' space is
        fprintf('[%0.3f,%0.3f] \n',color_markers(2,1),color_markers(2,2));
        
    else
        [color_markers, color_marker_labels] = histology_lab_colors('hem');
        ncolors = length(color_marker_labels);
    end
    
    
    % Step 3: Classify Each Pixel Using the Nearest Neighbor Rule
    % Each color marker now has an 'a*' and a 'b*' value. You can classify each pixel in the lab_fabric image by calculating the Euclidean distance between that pixel and each color marker. The smallest distance will tell you that the pixel most closely matches that color marker. For example, if the distance between a pixel and the red color marker is the smallest, then the pixel would be labeled as a red pixel.
    % Create an array that contains your color labels, i.e., 0 = background, 1 = red, 2 = green, 3 = purple, 4 = magenta, and 5 = yellow.
    
    color_labels = 0:ncolors-1;
    
    
    % Initialize matrices to be used in the nearest neighbor classification.
    a = double(a);
    b = double(b);
    distance = zeros([size(a), ncolors]);
    
    % Perform classification
    for count = 1:ncolors
        distance(:,:,count) = ( (a - color_markers(count,1)).^2 + (b - color_markers(count,2)).^2 ).^0.5;
    end
    
    [~, label] = min(distance,[],3);
    label = color_labels(label);
    clear distance;
    
    
    % Step 4: Display Results of Nearest Neighbor Classification
    % The label matrix contains a color label for each pixel in the fabric image. Use the label matrix to separate objects in the original fabric image by color.
    rgb_label = repmat(label,[1 1 3]);
    segmented_images = zeros([size(Im), ncolors],'uint8');
    
    for count = 1:ncolors
        color = Im;
        color(rgb_label ~= color_labels(count)) = 0;
        segmented_images(:,:,:,count) = color;
    end
    
    if false
        figure;
        subplot(1,4,1), imshow(segmented_images(:,:,:,1)), title('White objects');
        subplot(1,4,2), imshow(segmented_images(:,:,:,2)), title('Dark purple objects');
        subplot(1,4,3), imshow(segmented_images(:,:,:,3)), title('Light purple objects');
        subplot(1,4,4), imshow(segmented_images(:,:,:,4)), title('Transparent objects');
        subplot(1,4,4), imshow(segmented_images(:,:,:,5)), title('Red objects');
        
        figure;
        % subplot(1,3,1),
        imshow(segmented_images(:,:,:,3)+segmented_images(:,:,:,2)+segmented_images(:,:,:,4)), title('Light & purple objects');
        
        Im_NoHem = sum(segmented_images(:,:,:,1:4),4);
        
        imshowpair(Im,Im_NoHem,'montage')
        
        
    end
    
    
    % subplot(1,3,2), imfill(imclearborder(segmented_images(:,:,:,3)+segmented_images(:,:,:,2)),'hole'); title('Light & purple objects');
    %tissue_im = segmented_images(:,:,:,2)+segmented_images(:,:,:,3)+segmented_images(:,:,:,4);
    tissue_im = sum(segmented_images(:,:,:,2:end),4);
    
    tissue_im_filled = imfill(tissue_im,'holes');
    if false; imshowpair(Im,tissue_im_filled,'montage'); end
    
    
    %binary image of grayscale image
    imgb = imbinarize(rgb2gray(tissue_im));
    
    % Add back in dark nuclei
    if true
        dark_purple_part = segmented_images(:,:,:,2);
        dark_purple_part_gray = rgb2gray(dark_purple_part);
        dark_purple_part_bw = imbinarize(dark_purple_part_gray);
        imgb = imgb | dark_purple_part_bw;
    end
    
    % fill holes
    imgb_filled = imfill(~imgb,'holes');
    
    if false; imshowpair(Im,imgb_filled,'montage'); end
    
    % Dilate image to ty to minimize false holes
    se90 = strel('line', 3, 90);
    se0 = strel('line', 3, 0);
    imgb_filled_dilated = ~imdilate(~imgb_filled, [se90 se0]);
    if false; imshowpair(Im,imgb_filled_dilated,'montage'); end
    
    % Follow tutorial to fill small holes
    %https://blogs.mathworks.com/steve/2008/08/05/filling-small-holes/
    % Fill small holes
    
    % Step 1: Fill all holes using imfill:
    imtmp = ~imgb_filled_dilated;
    imtmp_filled = imfill(imtmp, 'holes'); % all holes filled
    if false; figure; imshow(imtmp_filled); end
    
    % Step 2: Identify the hole pixels using logical operators:
    holes = imtmp_filled & ~imtmp;
    % imshow(holes)
    title('Hole pixels identified')
    
    % Step 3: Use bwareaopen on the holes image to eliminate small holes:
    bigholes = bwareaopen(holes, 200);
    % imshow(~bigholes)
    title('Only the big holes')
    
    
    % Step 4: Use logical operators to identify small holes:
    smallholes = holes & ~bigholes;
    % imshow(smallholes)
    title('Only the small holes')
    
    % Step 5: Use a logical operator to fill in the small holes in the original image:
    Imholesfilled = ~(imtmp | smallholes);
    
    % Smooth the boundaries
    smoothing_factor =5; %disk radius
    ImSmoothed = imfilter(Imholesfilled,fspecial('disk',smoothing_factor));
    
    
    % Trim the edges
    cto = 50;
    ImTrimmed = Im(cto:end-cto,cto:end-cto,:);
    ImFinal = ImSmoothed(cto:end-cto,cto:end-cto,:);
    imshowpair(ImTrimmed,ImFinal,'montage')
    title('Original image vs initial condition')
    
    axis on; xticks([]);    yticks([])
    
else
    
    clc; clear all; close all;
    
    % Read the image file
    ImDir='e:/brandon/microscope_images/';
    ImFile = 'LTA1_x40_edge.jpg';%'LTA1_x40_edge_of_hem.jpg';
    Im = imread(sprintf('%s%s',ImDir,ImFile));
    
    % Convert RGB to L*a*b
    lab_im = rgb2lab(Im);
    
    %Calculate the mean 'a*' and 'b*' value for each area that you extracted with roipoly. These values serve as your color markers in 'a*b*' space.
    a = lab_im(:,:,2);
    b = lab_im(:,:,3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET COLORS %%%%%%%%%%%%%%%%%%%%
    
    % Define a function to grab a square for each color
    mysquare = @(x,y,s) ([x+floor(s/2).*[-1;-1;+1;+1;-1], y+floor(s/2).*[-1;+1;+1;-1;-1]]);
    
    % Rgioncoordinates for areas of each color (this is done manually right now)
    white_region = mysquare(762,503,14);
    dpurple_region = mysquare(849,565,10);
    ddpurple_region = mysquare(875,947,10);
    lpurple_region = mysquare(790,602,8);
    transparent_region = mysquare(1918,709,10);
    
    region_coordinates(:,:,1) = white_region;
    region_coordinates(:,:,2) = dpurple_region;
    region_coordinates(:,:,3) = lpurple_region;
    region_coordinates(:,:,4) = transparent_region;
    % region_coordinates(:,:,5) = ddpurple_region; %UNCOMMENTING THIS TO ACCOUNT FOR DARK PURPLE BREAKS EVERYTHING
    
    ncolors = 4;
    sample_regions = false([size(Im,1) size(Im,2) ncolors]);
    
    for count = 1:ncolors
        sample_regions(:,:,count) = roipoly(Im,region_coordinates(:,1,count), region_coordinates(:,2,count));
    end
    
    
    
    if false
        color_markers = zeros([ncolors, 2]);
        for count = 1:ncolors
            color_markers(count,1) = mean2(a(sample_regions(:,:,count)));
            color_markers(count,2) = mean2(b(sample_regions(:,:,count)));
        end
        
        %For example, the average color of the dark purple sample region in 'a*b*' space is
        fprintf('[%0.3f,%0.3f] \n',color_markers(2,1),color_markers(2,2));
        
    else
        [color_markers, color_marker_labels] = histology_lab_colors();
    end
    
    % Step 3: Classify Each Pixel Using the Nearest Neighbor Rule
    % Each color marker now has an 'a*' and a 'b*' value. You can classify each pixel in the lab_fabric image by calculating the Euclidean distance between that pixel and each color marker. The smallest distance will tell you that the pixel most closely matches that color marker. For example, if the distance between a pixel and the red color marker is the smallest, then the pixel would be labeled as a red pixel.
    % Create an array that contains your color labels, i.e., 0 = background, 1 = red, 2 = green, 3 = purple, 4 = magenta, and 5 = yellow.
    
    color_labels = 0:ncolors-1;
    
    
    % Initialize matrices to be used in the nearest neighbor classification.
    a = double(a);
    b = double(b);
    distance = zeros([size(a), ncolors]);
    
    % Perform classification
    for count = 1:ncolors
        distance(:,:,count) = ( (a - color_markers(count,1)).^2 + (b - color_markers(count,2)).^2 ).^0.5;
    end
    
    [~, label] = min(distance,[],3);
    label = color_labels(label);
    clear distance;
    
    
    % Step 4: Display Results of Nearest Neighbor Classification
    % The label matrix contains a color label for each pixel in the fabric image. Use the label matrix to separate objects in the original fabric image by color.
    rgb_label = repmat(label,[1 1 3]);
    segmented_images = zeros([size(Im), ncolors],'uint8');
    
    for count = 1:ncolors
        color = Im;
        color(rgb_label ~= color_labels(count)) = 0;
        segmented_images(:,:,:,count) = color;
    end
    
    if false
        figure;
        subplot(1,4,1), imshow(segmented_images(:,:,:,1)), title('White objects');
        subplot(1,4,2), imshow(segmented_images(:,:,:,2)), title('Dark purple objects');
        subplot(1,4,3), imshow(segmented_images(:,:,:,3)), title('Light purple objects');
        subplot(1,4,4), imshow(segmented_images(:,:,:,4)), title('Transparent objects');
        
        
        figure;
        % subplot(1,3,1),
        imshow(segmented_images(:,:,:,3)+segmented_images(:,:,:,2)+segmented_images(:,:,:,4)), title('Light & purple objects');
    end
    
    
    % subplot(1,3,2), imfill(imclearborder(segmented_images(:,:,:,3)+segmented_images(:,:,:,2)),'hole'); title('Light & purple objects');
    %tissue_im = segmented_images(:,:,:,2)+segmented_images(:,:,:,3)+segmented_images(:,:,:,4);
    tissue_im = sum(segmented_images(:,:,:,2:end),4);
    
    tissue_im_filled = imfill(tissue_im,'holes');
    if false; imshowpair(Im,tissue_im_filled,'montage'); end
    
    
    %binary image of grayscale image
    imgb = imbinarize(rgb2gray(tissue_im));
    
    % Add back in dark nuclei
    if true
        dark_purple_part = segmented_images(:,:,:,2);
        dark_purple_part_gray = rgb2gray(dark_purple_part);
        dark_purple_part_bw = imbinarize(dark_purple_part_gray);
        imgb = imgb | dark_purple_part_bw;
    end
    
    % fill holes
    imgb_filled = imfill(~imgb,'holes');
    
    if false; imshowpair(Im,imgb_filled,'montage'); end
    
    % Dilate image to ty to minimize false holes
    se90 = strel('line', 3, 90);
    se0 = strel('line', 3, 0);
    imgb_filled_dilated = ~imdilate(~imgb_filled, [se90 se0]);
    if false; imshowpair(Im,imgb_filled_dilated,'montage'); end
    
    % Follow tutorial to fill small holes
    %https://blogs.mathworks.com/steve/2008/08/05/filling-small-holes/
    % Fill small holes
    
    % Step 1: Fill all holes using imfill:
    imtmp = ~imgb_filled_dilated;
    imtmp_filled = imfill(imtmp, 'holes'); % all holes filled
    if false; figure; imshow(imtmp_filled); end
    
    % Step 2: Identify the hole pixels using logical operators:
    holes = imtmp_filled & ~imtmp;
    % imshow(holes)
    title('Hole pixels identified')
    
    % Step 3: Use bwareaopen on the holes image to eliminate small holes:
    bigholes = bwareaopen(holes, 200);
    % imshow(~bigholes)
    title('Only the big holes')
    
    
    % Step 4: Use logical operators to identify small holes:
    smallholes = holes & ~bigholes;
    % imshow(smallholes)
    title('Only the small holes')
    
    % Step 5: Use a logical operator to fill in the small holes in the original image:
    Imholesfilled = ~(imtmp | smallholes);
    
    % Smooth the boundaries
    smoothing_factor =5; %disk radius
    ImSmoothed = imfilter(Imholesfilled,fspecial('disk',smoothing_factor));
    
    
    % Trim the edges
    cto = 50;
    ImTrimmed = Im(cto:end-cto,cto:end-cto,:);
    ImFinal = ImSmoothed(cto:end-cto,cto:end-cto,:);
    imshowpair(ImTrimmed,ImFinal,'montage')
    title('Original image vs initial condition')
    
    axis on; xticks([]);    yticks([])
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if false
    
    %%
    % clear the workspace
    clc; clear all; close all;
    
    % Read the image file
    ImDir='e:/brandon/microscope_images/';
    ImFile = 'LTA1_x40_edge.jpg';
    Im = imread(sprintf('%s%s',ImDir,ImFile));
    
    % Thresholding the image on each color plane
    Im = im2double(Im);
    [r, c, p] = size(Im);
    
    imR = squeeze(Im(:,:,1));
    imG = squeeze(Im(:,:,2));
    imB = squeeze(Im(:,:,3));
    
    %Thresholding on individual planes
    if false
        imBinaryR = im2bw(imR, graythresh(imR));
        imBinaryG = im2bw(imG, graythresh(imG));
        imBinaryB = im2bw(imB, graythresh(imB));
    else
        imBinaryR = imbinarize(imR);
        imBinaryG = imbinarize(imG);
        imBinaryB = imbinarize(imB);
    end
    
    imBinary = imcomplement(imBinaryR&imBinaryG&imBinaryB);
    
    imshow(imBinary)
    
    % Morphological opening (may not work hear)
    se = strel('disk', 7);
    imCleanDisks = imopen(imBinary, se);
    
    % Fill holes and clear boarder
    imFilled = imfill(imBinary, 'holes');
    imClear = imclearborder(imBinary);
    
    % Count the number of objects in an image
    [labels, numLabels] = bwlabel(imCleanDisks);
    disp(['Number of objects detected: ' num2str(numLabels)]);
    
    % Initizlize matrices
    rLabel = zeros(r,c);
    gLabel = zeros(r,c);
    bLabel = zeros(r,c);
    
    % Get average color vector for each labeled region
    for ii = 1:numLabels
        rLabel(labels==ii) = median(imR(labels==ii));
        gLabel(labels==ii) = median(imG(labels==ii));
        bLabel(labels==ii) = median(imB(labels==ii));
    end
    imLabel = cat(3, rLabel, gLabel, bLabel);
    
    imshow(imLabel);
    impixelinfo(gcf);
    
    % get desired color to segment
    [x, y] = ginput(1);
    selColor = imLabel(floor(y), floor(x), :);
    
    % Convert to LAB color space
    C = makecform('srgb2lab');
    imLAB = applycform(imLabel,C);
    imSelLab=applycform(selColor,C);
    
    % Extract a* and b* values
    imA = imLAB(:,:,2);
    imB = imLAB(:,:,3);
    imSelA = imSelLab(1,2); %extract a*
    imSelB = imSelLab(1,3); %extract a*
    
    distThresh = 10;
    imMask = zeros(r,c);
    imDist = hypot(imA-imSelA,imB-imSelB);
    imMask(imDist<distThresh) = 1;
    [cLabel, cNum] = bwlabel(imMask);
    imSeg = repmat(selColor, [r,c,1]).*repmat(imMask, [1, 1, 3]);
    imshow(imSeg)
    
    
    
    % Plot the images
    figure
    subplot(2,2, 1)
    imshow(Im)
    title('Original Image')
    
    subplot(2,2, 2)
    imshow(imBinary)
    title('Binary Image')
    
    subplot(2,2, 3)
    imshow(imFilled)
    title('Filled Binary Image')
    
    subplot(2,2, 4)
    imshow(imCleanDisks)
    title('Disks Image')
    
    
    
    
    %%
    
    se = strel('disk',1);
    tissue_im_closed = imclose(tissue_im,se);
    imshow(tissue_im_closed)
    
    if false
        % disk average filter
        % myfilter = fspecial('disk',3);
        % tissue_im_filtered = imfilter(tissue_im,myfilter);
    elseif false
        %median filter
        % tissue_im_filtered = medfilt2(rgb2gray(tissue_im));
        figure
        imshow(tissue_im_filtered)
    end
    
    % use super pixels
    [L,N] = superpixels(tissue_im,500);
    figure
    BW = boundarymask(L);
    imshow(imoverlay(tissue_im,BW,'cyan'),'InitialMagnification',67)
    
    % Further process the white parts
    figure
    whitefilled = imfill(segmented_images(:,:,:,1),'holes');
    
    
    % The binary gradient mask shows lines of high contrast in the image. These lines do not quite delineate the outline of the object of interest. Compared to the original image, you can see gaps in the lines surrounding the object in the gradient mask. These linear gaps will disappear if the Sobel image is dilated using linear structuring elements, which we can create with the strel function.
    se90 = strel('line', 3, 90);
    se0 = strel('line', 3, 0);
    
    % The binary gradient mask is dilated using the vertical structuring element followed by the horizontal structuring element. The imdilate function dilates the image.
    whitefilled_dil = imdilate(whitefilled, [se90 se0]);
    imshow(whitefilled_dil)
    %% Use super pixesl
    
    clc; clear all; close all;
    
    % Read the image file
    ImDir='e:/brandon/microscope_images/';
    ImFile = 'LTA1_x40_edge.jpg';
    Im = imread(sprintf('%s%s',ImDir,ImFile));
    % imshow(Im)
    
    
    % Calculate superpixels of the image.
    % [L,N] = superpixels(Im,1000);
    [L,N] = superpixels(segmented_images(:,:,:,3)+segmented_images(:,:,:,2),1000);
    
    % Display the superpixel boundaries overlaid on the original image.
    figure
    BW = boundarymask(L);
    imshow(imoverlay(Im,BW,'cyan'),'InitialMagnification',67)
    
    %%
    clc; clear all; close all;
    
    % Read the image file
    ImDir='e:/brandon/microscope_images/';
    ImFile = 'LTA1_x40_edge.jpg';
    Im = imread(sprintf('%s%s',ImDir,ImFile));
    
    Im = im2double(Im);
    Img = rgb2gray(Im);
    
    [BWs, threshold] = edge(Img, 'sobel');
    % figure, imshow(BWs), title('binary gradient mask');
    fudgeFactor = 1.0;
    BWs = edge(Img,'sobel', threshold * fudgeFactor);
    % figure, imshow(BWs), title('binary gradient mask');
    
    % Step 3: Dilate the Image
    % The binary gradient mask shows lines of high contrast in the image. These lines do not quite delineate the outline of the object of interest. Compared to the original image, you can see gaps in the lines surrounding the object in the gradient mask. These linear gaps will disappear if the Sobel image is dilated using linear structuring elements, which we can create with the strel function.
    se90 = strel('line', 5, 90);
    se0 = strel('line', 5, 0);
    
    % The binary gradient mask is dilated using the vertical structuring element followed by the horizontal structuring element. The imdilate function dilates the image.
    BWsdil = imdilate(BWs, [se90 se0]);
    % figure, imshow(BWsdil), title('dilated gradient mask');
    
    
    % Step 4: Fill Interior Gaps
    % The dilated gradient mask shows the outline of the cell quite nicely, but there are still holes in the interior of the cell. To fill these holes we use the imfill function.
    BWdfill = imfill(BWsdil, 'holes');
    % figure, imshow(BWdfill);
    % title('binary image with filled holes');
    
    % Step 5: Remove Connected Objects on Border
    % The cell of interest has been successfully segmented, but it is not the only object that has been found. Any objects that are connected to the border of the image can be removed using the imclearborder function. The connectivity in the imclearborder function was set to 4 to remove diagonal connections.
    % BWnobord = imclearborder(BWdfill);
    % figure, imshow(BWnobord), title('cleared border image');
    
    % Step 6: Smoothen the Object
    % Finally, in order to make the segmented object look natural, we smoothen the object by eroding the image twice with a diamond structuring element. We create the diamond structuring element using the strel function.
    seD = strel('diamond',1);
    BWnobord = BWdfill;
    BWfinal = imerode(BWnobord,seD);
    BWfinal = imerode(BWfinal,seD);
    figure, imshow(BWfinal), title('segmented image');
    
    % figure;
    nsubfig = 2;
    subplot(1,nsubfig,1), imshow(Im); title('Original Image');% original
    % subplot(1,nsubfig,2), imshow(Img); title('Grayscale Image');% gray image
    % subplot(1,nsubfig,3), imshow(~BWs); title('Binary gradient mask');
    % subplot(1,nsubfig,4), imshow(~BWsdil); title('Dilated gradient mask');
    subplot(1,nsubfig,2), imshow(~BWdfill); title('binary image with filled holes');
    % subplot(1,nsubfig,6), imshow(~BWfinal), title('segmented image');
    
    
    %%
    clc; clear all; close all;
    
    % Read the image file
    ImDir='e:/brandon/microscope_images/';
    ImFile = 'LTA1_x40_edge.jpg';
    Im = imread(sprintf('%s%s',ImDir,ImFile));
    
    
    % spacial average of image
    filter2(fspecial('average',3),J)/255;
    
    %% Remove objects smaller than threshold number of pixels
    
    clc; clear all; close all;
    
    % Read the image file
    ImDir='e:/brandon/microscope_images/';
    ImFile = 'LTA1_x40_edge.jpg';
    Im = imread(sprintf('%s%s',ImDir,ImFile));
    
    
    % spacial average of image
    Img = rgb2gray(Im);
    
    BW = imbinarize(Img);
    
    
    % Remove small artifacts (fewer than P pixels)
    P = 1000;
    BW2 = bwareaopen(BW,P,8);
    
    figure;
    imshowpair(BW,BW2)
    % subplot(1,2,1), imshow(BW)
    % subplot(1,2,2), imshow(BW2)
    
    
    %% morphologically close picure
    clc; clear all; close all;
    
    % Read the image file
    ImDir='e:/brandon/microscope_images/';
    ImFile = 'LTA1_x40_edge.jpg';
    Im = imread(sprintf('%s%s',ImDir,ImFile));
    
    
    % spacial average of image
    Img = rgb2gray(Im);
    
    % BW = imbinarize(Img);
    BW = imbinarize(Img,'adaptive','Sensitivity',0.95);
    
    se = strel('disk',5);
    
    ImClosed = imclose(~BW,se);
    % imshowpair(Im,ImClosed,'montage')
    
    ImFillClosed= imfill(ImClosed,'holes');
    
    imshowpair(Im,~ImFillClosed,'montage')
    
end

