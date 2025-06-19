clear all;
close all;
clc;


%ex. 1
% basically go over the whole image with that mask
thetas = [30,45,60];
filters = zeros(3, 3, 4);
filters(:,:,1) = [-1 -1 -1; 2 2 2; -1 -1 -1]; % horiz
filters(:,:,2) = [-1 -1 2; -1 2 -1; 2 -1 -1]; % 45 deg
filters(:,:,3) = [-1 2 -1; -1 2 -1; -1 2 -1]; % vert
filters(:,:,4) = [2 -1 -1; -1 2 -1; -1 -1 2]; % -45 deg
DisplayLines(thetas, 1, filters);
DisplayLines(thetas, 10, filters);

%ex2
close all;
clear all;
clc;

image = imread("SourceImages\dodecahedron_example.png");
minT = min(image(:));
maxT = max(image(:));
T = (minT + maxT) / 2;
threshhold = 0.001;
converged = false;
while not(converged)
    smlrThanT = image(image < T);
    bigThanT = image(image>T);
    meanSmaller = mean(smlrThanT);
    meanBigger = mean(bigThanT);
    newT = (meanSmaller + meanBigger) / 2;
    if abs(T-newT) < threshhold
       converged = true;
    end
    T = newT;
end
T
image = (image > T);
imshow(image);

%ex. 3
clc;
close all;
clear all;
image = im2double(imread("SourceImages\medtest.png"));
%f1 = figure;
% imshow(image, []);
% impixelinfo;
pixelsOfInterest = [
    [213, 251]; %vertibrade?
    [190, 253]; % ribs? % rib1
    [127, 255]; % rib2
    [132, 252];
    [101, 212]; %rib3
    [105, 209];

    [102, 164]; % rib4
    [99, 160];

    [108, 117];
    [112, 120];

    [155, 79];
    [144, 91];

    [240, 79];

    [353, 111];

    [368, 157];
    [373, 152];

    [364, 201];
    [369, 201];

    [349, 243];
    [351, 246];

    [358, 109]; % extra seed on dim rib

    [296, 271];
    [262, 253];

    [309, 248]; %shroom cancer
];
threshhold = 0.08;
rgbImage = repmat(image, [1 1 3]);

for seedIndex = 1:size(pixelsOfInterest, 1)
    y = pixelsOfInterest(seedIndex, 1);
    x = pixelsOfInterest(seedIndex, 2);

    region = regiongrowing(image, x, y, threshhold);

    redChannel = rgbImage(:,:,1);
    greenChannel = rgbImage(:,:,2);
    blueChannel = rgbImage(:,:,3);

    redChannel(region) = 1;     
    greenChannel(region) = 0;    
    blueChannel(region) = 0;     

    rgbImage(:,:,1) = redChannel;
    rgbImage(:,:,2) = greenChannel;
    rgbImage(:,:,3) = blueChannel;
end

figure;
imshow(rgbImage);
title("Region Growing Overlay in Red");
impixelinfo;

% 3b
f2 = figure;
imhist(image);
title("Histogram of the whole image");

figure;
snippet = image(200:300, 180:280);
snippet = snippet > 0.9;
imshow(snippet)
title("Threshholded (>0.9) snippet of image")

function DisplayLines(thetas, lineSize, filters)
    filterNames = [
      "Horizontal", "45 degrees", "Vertical", "-45 degrees"
    ];
    A = zeros(500); [j, i] = find(A<Inf);
    for index = 1:length(thetas) 
        q = find( abs(j-tan(thetas(index)*pi/180)*i)<=lineSize);
        A(q) = 100;
    end    
    fig = figure;
    for i = 1:4
        filtered = imfilter(A, filters(:,:,i));
        subplot(2, 2, i);
        imshow(filtered);
        title(sprintf('Filter %s - %dpx', filterNames(i), lineSize));
    end
end



