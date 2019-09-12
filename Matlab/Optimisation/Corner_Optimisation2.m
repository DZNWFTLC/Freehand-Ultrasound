function output = Corner_Optimisation2(T)
global Window;
Mark = rgb2gray(imread('4.png'));
% figure(2)
% imshow(Window);
tform = projective2d(T);
%Window = imbilatfilt(Window);
Template = imwarp(imcomplement(Mark),tform, 'SmoothEdges', true);
Col_valid = any(Template,1);
Row_valid = any(Template,2);
Template = imcomplement(Template(Row_valid,Col_valid));
sizeT = size(Template);
sizeW = size(Window);
scale = min((sizeW-20)./sizeT);
Template = imresize(Template,scale);


c = normxcorr2(Template,Window);
output = 1- max(c(:));
max(c(:))
