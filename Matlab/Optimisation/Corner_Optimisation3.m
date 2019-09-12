function output = Corner_Optimisation3(offset)
global Window;
global Coord;
Mark = rgb2gray(imread('4.png'));
% figure(2)
% imshow(Window);
[Row,Col] = size(Mark);
Coord_ori = [0 0; Row 0; Row Col;0 Col;];
offset_index = reshape(offset,[4,2]);
Coord_new = Coord + offset_index;
tform = estimateGeometricTransform(Coord_ori,Coord_new,'projective');
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
output = max(c(:));
max(c(:))
