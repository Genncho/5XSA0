%  ==========================================================================
%  | Instructions to test individual exercise:                              |
%  | Add a return at the end of the exercise code, so it stops executing    |
%  |    further code                                                        |
%  ==========================================================================

% ex1
clear all;
close all;
clc;
image = zeros(200, 200, "uint8");
x_start = 50;
y_start = 50;
width = 100;
height = 100;

image(y_start:y_start+height-1, x_start:x_start+width-1) = 255;
figure; imshow(image, []); title("White rectangle on a black background");

fspec = fft2(image);
figure; imshow(abs(fspec), []); title("Fourier spectrum of the rectangle");

fspec_shifted = fftshift(fspec);
figure; imshow(abs(fspec_shifted), []); title("Centered Fourier spectrum");

magnitude_spectrum = log(1 + abs(fspec_shifted));
figure; imshow(abs(magnitude_spectrum), []); title("Visually enchanced spectrum");

% ex2
clear all;
close all;
clc;

snapshot = imread("files\mri_snapshot.jpg");
snapshot = rgb2gray(snapshot);
figure; imshow(snapshot, []); title("Initial image");

[M,N] = size(snapshot);
freqSpec = fftshift(fft2(snapshot));
[u, v] = meshgrid(-floor(N/2):ceil(N/2)-1, -floor(M/2):ceil(M/2)-1);
D = sqrt(u.^2 + v.^2);
D_max = max(D(:));
cutoff_percent = 0.1;
D0 = cutoff_percent * D_max;
H = double(D <= D0);
figure; imshow(H, []); title("Filtered Transfer function");

figure; mesh(u, v, H); title("The Transfer function as 3D mesh");
G = H .* freqSpec;
img_filtered = real(ifft2(ifftshift(G)));
figure; imshow(uint8(img_filtered)); title("Resulting image from filtering");

% ex3
clear all;
close all;
clc;

snapshot = imread("files\mri_snapshot.jpg");
snapshot = rgb2gray(snapshot);
figure; imshow(snapshot, []);  title("Initial image")

[M,N] = size(snapshot);
freqSpec = fftshift(fft2(snapshot));
[u, v] = meshgrid(-floor(N/2):ceil(N/2)-1, -floor(M/2):ceil(M/2)-1);
D = sqrt(u.^2 + v.^2);
D_max = max(D(:));
cutoff_percent = 0.1;
D0 = cutoff_percent * D_max;
H = 1 - double(D <= D0);
figure; imshow(H, []); title("Filtered Transfer function");

figure; mesh(u, v, H); title("The Transfer function as 3D mesh");

G = H .* freqSpec;
img_filtered = real(ifft2(ifftshift(G)));
figure; imshow(img_filtered, []); title("Resulting image from filtering");

%ex 4 
clear all;
close all;
clc;

snapshot = imread("files\mri_snapshot.jpg");
snapshot = rgb2gray(snapshot);

[M, N] = size(snapshot);
freqSpec = fftshift(fft2(snapshot));
[u, v] = meshgrid(-floor(N/2):ceil(N/2)-1, -floor(M/2):ceil(M/2)-1);
D = sqrt(u.^2 + v.^2);
D_max = max(D(:));

cutoff_percents = [0.01, 0.05, 0.5];

figure;
for i = 1:length(cutoff_percents)
    cutoff_percent = cutoff_percents(i);
    D0 = cutoff_percent * D_max;
    H = double(D <= D0);

    G = H .* freqSpec;
    img_filtered = real(ifft2(ifftshift(G)));

    subplot(3, 3, (i-1)*3 + 1);
    imshow(H, []);
    title(sprintf('H (D0=%.2f)', cutoff_percent));

    subplot(3, 3, (i-1)*3 + 2);
    mesh(u, v, H);
    title(sprintf('Mesh (D0=%.2f)', cutoff_percent));
    axis tight; view(30, 30);

    subplot(3, 3, (i-1)*3 + 3);
    imshow(img_filtered, []);
    title(sprintf('Filtered (D0=%.2f)', cutoff_percent));
end

return
%ex 5 
clear all;
close all;
clc;

snapshot = imread("files\mri_snapshot.jpg");
snapshot = rgb2gray(snapshot);
[M, N] = size(snapshot);
freqSpec = fftshift(fft2(snapshot));
[u, v] = meshgrid(-floor(N/2):ceil(N/2)-1, -floor(M/2):ceil(M/2)-1);
D = sqrt(u.^2 + v.^2);
D_max = max(D(:));

cutoff_vals = [0.01, 0.05, 0.5];
figure('Name','Task 5: Highpass Filtering','NumberTitle','off');

for i = 1:length(cutoff_vals)
    D0 = cutoff_vals(i) * D_max;
    H = 1 - double(D <= D0);  
    G = H .* freqSpec;       
    img_filtered = real(ifft2(ifftshift(G)));  

    subplot(3, 3, (i-1)*3 + 1);
    imshow(H, []);
    title(sprintf('HPF H, D0 = %.2f', cutoff_vals(i)));

    subplot(3, 3, (i-1)*3 + 2);
    mesh(u, v, H);
    title(sprintf('Mesh HPF, D0 = %.2f', cutoff_vals(i)));
    axis tight; view(30, 30);

    subplot(3, 3, (i-1)*3 + 3);
    imshow(img_filtered, []);
    title(sprintf('HPF Img, D0 = %.2f', cutoff_vals(i)));
end



