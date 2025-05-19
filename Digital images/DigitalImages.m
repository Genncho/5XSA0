%  ==========================================================================
%  | Instructions to test individual exercise:                              |
%  | Add a return at the end of the exercise code, so it stops executing    |
%  |    further code                                                        |
%  ==========================================================================


% ex1
clear all;
close all;
clc;
image = imread("files\chestxray.tif");
imshow(image, []);
title('Chest Xray image')
impixelinfo;
imwrite(image, "savedFiles\chestXraySaved.tif");


% ex2
clear all;
close all;
clc;
m1 = [1 2 3; 4 5 6; 7 8 9];
m2 = [10 1 4; 4 5 2; 7 2 3];
m3 = m1 * m2;
m4 = m1 .* m2;

disp('Matrix m1:');
disp(m1);

disp('Matrix m2 (Identity matrix):');
disp(m2);

disp('Matrix m3 = m1 * m2 (Matrix multiplication):');
disp(m3);

disp('Matrix m4 = m1 .* m2 (Element-wise multiplication):');
disp(m4);

% ex3
clear all;
close all;
clc;

cropPx = 50;
image = imread("files\chestxray.tif");
flippedImage = image(:, end:-1:1, :);
croppedImage = flippedImage(cropPx:end-cropPx, cropPx:end-cropPx);
subsampledImage = croppedImage(1:2:end, 1:2:end);

figure; imshow(flippedImage, []); title('Flipped Horizontally');
figure; imshow(croppedImage, []); title('Flipped + Cropped');
figure; imshow(subsampledImage, []); title('Flipped + Cropped + Subsampled');

% ex4
clear all;
close all;
clc;
image = imread("files\mri_snapshot.jpg");
grayImage = rgb2gray(image);
equalizedImage = histeq(grayImage);

avgFilter3x3 = fspecial("average", [3 3]);
avgFilter7x7 = fspecial("average", [7 7]);

avgFImage3 = imfilter(equalizedImage, avgFilter3x3);
avgFImage7 = imfilter(equalizedImage, avgFilter7x7);

medFilteredImage3 = medfilt2(equalizedImage, [3 3]);
medFilteredImage7 = medfilt2(equalizedImage, [7 7]);

figure; imshow(grayImage); title("Grayscale Image");
figure; imhist(image); title("Grayscale Image Histogram");
figure; imshow(equalizedImage); title("Image After Histogram Equalization")
figure; imhist(equalizedImage); title('Histogram After Equalization');
figure; imshow(avgFImage3); title('Equalized image after 3x3 averaging filter');
figure; imshow(avgFImage7); title('Equalized image after 7x7 averaging filter');
figure; imshow(medFilteredImage3); title('Equalized image after 3x3 median filter');
figure; imshow(medFilteredImage7); title('Equalized image after 7x7 median filter');


return;
% extra code for analyzing differentces between the filters
subtractedImage = imsubtract(avgFImage3, medFilteredImage3);

% Display the result
imshow(subtractedImage, []);
title('Average (3x3) - Median (3x3) Filtered Image');




