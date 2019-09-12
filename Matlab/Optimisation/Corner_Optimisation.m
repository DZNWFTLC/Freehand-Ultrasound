function output = Corner_Optimisation(Corners_mat,CameraImage)
Mark = rgb2gray(imread('4.png'));
Coord_Col = Corners_mat(1:4);
Coord_Row = Corners_mat(5:end);
Range = [min(Coord_Row),max(Coord_Row);
    min(Coord_Col),max(Coord_Col)];
Window = CameraImage(Range(1,1)-10:Range(1,2)+10,Range(2,1)-10:Range(2,2)+10);
% figure(2)
% imshow(Window);

%Window = imbilatfilt(Window);
Coord = [Coord_Col.'-min(Coord_Col) Coord_Row.'-min(Coord_Row)];
[Row,Col] = size(Mark);
Coord_ori = [0 0; Row 0; Row Col;0 Col;];
tform = estimateGeometricTransform(Coord_ori,Coord,'projective');
Template = imwarp(imcomplement(Mark),tform, 'SmoothEdges', true);
Col_valid = any(Template,1);
Row_valid = any(Template,2);
Template = imcomplement(Template(Row_valid,Col_valid));
% figure(3)
% imshow(Template);
%     ratio_x =  / rect.width;
% 	ratio_y = (float)dist_y / rect.height;
%
% 	if (abs(ratio_x - 1) <= abs(ratio_y - 1))
% 		Resize = Size(round(ratio_x * rect.width), round(ratio_x * rect.height));
% 	else
% 		Resize = Size(round(ratio_y * rect.width), round(ratio_y * rect.height));
% 	resize(img_perspective(rect), img_marker, Resize, 0, 0, INTER_AREA);

%imresize

c = normxcorr2(Template,Window);

figure, surf(c), shading flat
[ypeak, xpeak] = find(c==max(c(:)));

yoffSet = ypeak-size(Template,1);
xoffSet = xpeak-size(Template,2);

yoffSet = gather(ypeak-size(Template,1));
xoffSet = gather(xpeak-size(Template,2));
output = max(c(:));
figure
imshow(Window);
imrect(gca, [xoffSet+1, yoffSet+1, size(Template,2), size(Template,1)]);