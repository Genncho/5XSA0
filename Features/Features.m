%  ==========================================================================
%  | Instructions to test individual exercise:                              |
%  | Add a return at the end of the exercise code, so it stops executing    |
%  |    further code                                                        |
%  ==========================================================================

% ex. 1
clc;
clear all;
close all;
image = im2double(imread("InputImages\lena.bmp"));

% noise - for analysis
h = fspecial('gaussian', [5 5], 1); 
image_smooth = imfilter(image, h, 'replicate');
noise_estimate = image - image_smooth;

% Roberts
Gx_r = [1 0; 0 -1];
Gy_r = [0 1; -1 0];
convX_r = imfilter(image, Gx_r, 'replicate');
convY_r = imfilter(image, Gy_r, 'replicate');
roberts_result = sqrt(convX_r.^2 + convY_r.^2);

% Prewitt
Gx_p = [1 0 -1; 1 0 -1; 1 0 -1];
Gy_p = [1 1 1; 0 0 0; -1 -1 -1];
convX_p = imfilter(image, Gx_p, 'replicate');
convY_p = imfilter(image, Gy_p, 'replicate');
prewitt_result = sqrt(convX_p.^2 + convY_p.^2);

% Sobel
Gx_s = [1 0 -1; 2 0 -2; 1 0 -1];
Gy_s = [1 2 1; 0 0 0; -1 -2 -1];
convX_s = imfilter(image, Gx_s, 'replicate');
convY_s = imfilter(image, Gy_s, 'replicate');
sobel_result = sqrt(convX_s.^2 + convY_s.^2);

roberts_result = mat2gray(roberts_result);
prewitt_result = mat2gray(prewitt_result);
sobel_result = mat2gray(sobel_result);

figure;

subplot(2,3,1); imshow(image, []); title('Original Image');
subplot(2,3,2); imshow(noise_estimate, []); title('Noise');
subplot(2,3,4); imshow(roberts_result, []); title('Roberts Edge Detection');
subplot(2,3,5); imshow(prewitt_result, []); title('Prewitt Edge Detection');
subplot(2,3,6); imshow(sobel_result, []); title('Sobel Edge Detection');

figure
diff_sp = abs(sobel_result - prewitt_result);
%diff_sp = mat2gray(diff_sp);
imshow(diff_sp);
title('Sobel - Prewitt difference - direct subtraction result');

% noise is around the edges, so Sobel's noise reduction is barely noticable
% compared to Prewitt result

%ex. 2
clc;
clear all;
close all;

img_size = 200;
image = im2double(zeros(img_size, img_size, "uint8"));
center_x = 100;
center_y = 100;
radius = 50;
[x, y] = meshgrid(1:img_size, 1:img_size);
circle_mask = (x - center_x).^2 + (y - center_y).^2 <= radius^2;
image(circle_mask) = 255;
figure; imshow(image, []); title("White rectangle on a black background");

thresholds = [30/255, 70/255, 130/255];

% Roberts
Gx = [ 1 0 ; 0 -1 ];
Gy = [ 0 1 ; -1 0 ];

convX = imfilter(image, Gx, 'replicate');
convY = imfilter(image, Gy, 'replicate');
roberts_mag = sqrt(convX.^2 + convY.^2);

% Prewitt
Gx = [1 0 -1; 1 0 -1; 1 0 -1];
Gy = [1 1 1; 0 0 0; -1 -1 -1];
convX = imfilter(image, Gx, 'replicate');
convY = imfilter(image, Gy, 'replicate');
prewitt_mag = sqrt(convX.^2 + convY.^2);

figure;
for i = 1:length(thresholds)
    t = thresholds(i);
    rob_bin = roberts_mag > t;
    pre_bin = prewitt_mag > t;
    subplot(length(thresholds), 3, (i-1)*3 + 1);
    imshow(image); title('Original');
    subplot(length(thresholds), 3, (i-1)*3 + 2);
    imshow(rob_bin); title(['Roberts Threshold ', num2str(t*255)]);
    subplot(length(thresholds), 3, (i-1)*3 + 3);
    imshow(pre_bin); title(['Prewitt Threshold ', num2str(t*255)]);
end
hold on;

% The generated shape is too clean to show any meaningful difference in
% noise reduction or threshholding by applying Roberts or Prewitt filters.
% However with Lena, it becomes much more apperant the effect of the
% techniques and their threshholding values
clc;
clear all;
close all;

% ex 2.5
% with Lena as image
image = im2double(imread("InputImages\lena.bmp"));
thresholds = [30/255, 70/255, 130/255];

% Roberts
Gx = [ 1 0 ; 0 -1 ];
Gy = [ 0 1 ; -1 0 ];

convX = imfilter(image, Gx, 'replicate');
convY = imfilter(image, Gy, 'replicate');
roberts_mag = sqrt(convX.^2 + convY.^2);

% Prewitt
Gx = [1 0 -1; 1 0 -1; 1 0 -1];
Gy = [1 1 1; 0 0 0; -1 -1 -1];
convX = imfilter(image, Gx, 'replicate');
convY = imfilter(image, Gy, 'replicate');
prewitt_mag = sqrt(convX.^2 + convY.^2);

figure;
for i = 1:length(thresholds)
    t = thresholds(i);

    % Apply binary threshold
    rob_bin = roberts_mag > t;
    pre_bin = prewitt_mag > t;

    % Display images
    subplot(length(thresholds), 3, (i-1)*3 + 1);
    imshow(image); title('Original');

    subplot(length(thresholds), 3, (i-1)*3 + 2);
    imshow(rob_bin); title(['Roberts Threshold ', num2str(t*255)]);

    subplot(length(thresholds), 3, (i-1)*3 + 3);
    imshow(pre_bin); title(['Prewitt Threshold ', num2str(t*255)]);
end
hold on;

% ex. 3
clc;
clear all;
close all;

image_rgb = imread("InputImages\endoscopic_image.jpg");
image = im2double(rgb2gray(image_rgb));
threshhold = 50/255;

% Roberts
Gx_r = [1 0; 0 -1];
Gy_r = [0 1; -1 0];
convX_r = imfilter(image, Gx_r, 'replicate');
convY_r = imfilter(image, Gy_r, 'replicate');
roberts_result = sqrt(convX_r.^2 + convY_r.^2);
roberts_edges = roberts_result > threshhold;

% Prewitt
Gx_p = [1 0 -1; 1 0 -1; 1 0 -1];
Gy_p = [1 1 1; 0 0 0; -1 -1 -1];
convX_p = imfilter(image, Gx_p);
convY_p = imfilter(image, Gy_p);
prewitt_result = sqrt(convX_p.^2 + convY_p.^2);
prewitt_edges = prewitt_result > threshhold;

% Sobel
Gx_s = [1 0 -1; 2 0 -2; 1 0 -1];
Gy_s = [1 2 1; 0 0 0; -1 -2 -1];
convX_s = imfilter(image, Gx_s);
convY_s = imfilter(image, Gy_s);
sobel_result = sqrt(convX_s.^2 + convY_s.^2);
sobel_edges = sobel_result > threshhold;

% Canny
canny_edges = edge(image, "canny", [0.1 0.125]);

figure; imshow(image_rgb); title('Original RGB Image');
figure; imshow(image); title('Grayscale Image');
figure; imshow(roberts_edges); title('Roberts Edge Detection');
figure; imshow(prewitt_edges); title('Prewitt Edge Detection');
figure; imshow(sobel_edges); title('Sobel Edge Detection');
figure; imshow(canny_edges); title('Canny Edge Detection');

f2 =figure;
imshow(image_rgb); hold on;
visboundaries(canny_edges, 'Color', 'b', 'LineWidth', 0.7);
title('Canny Edges Overlayed on Original');

exportgraphics(f2,"outline.png")























