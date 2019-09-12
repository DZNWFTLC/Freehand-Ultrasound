%image 
%Transform matrix
%threshold 
clear
threshold = 200;
image = rgb2gray(imread('Kidney0900.bmp'));
pixel_valid = zeros(size(image));
[row,col] = find(image > threshold);
value = diag(image(row,col));
coordinate = [row col zeros(length(row),1)];