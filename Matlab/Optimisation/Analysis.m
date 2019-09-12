% mha1 = mha_read_header('Volume_test1.mha');
% mha2 = mha_read_header('Volume_test2.mha');
% 
% PC1 = mha_read_volume(mha1);
% PC2 = mha_read_volume(mha2);
% PC1 = Input_matrix(PC1,110,false);
% PC2 = Input_matrix(PC2,110,false);
% 
% figure(1)
% [Face,Vertices] = Marching_Cubes_Poisson(PC1,1);
% trisurf(Face,Vertices(:,1),Vertices(:,2),Vertices(:,3))
% axis equal
% grid on
% figure(2)
% [Face,Vertices] = Marching_Cubes_Poisson(PC2,1);
% trisurf(Face,Vertices(:,1),Vertices(:,2),Vertices(:,3))
% axis equal
% grid on

% mha1 = mha_read_header('Volume3_2.mha');
mha2 = mha_read_header('..\1_6.mha');

% PC1 = mha_read_volume(mha1);
PC2 = mha_read_volume(mha2);
% PC1 = Input_matrix(PC1,130,false);
PC2 = Input_matrix(PC2,500,false);

hold on
% [Face,Vertices] = Marching_Cubes_Poisson(PC1,1);
% trisurf(Face,Vertices(:,1),Vertices(:,2),Vertices(:,3))
axis equal
grid on
[Face,Vertices] = Marching_Cubes_Poisson(PC2,2);
trisurf(Face,Vertices(:,1),Vertices(:,2),Vertices(:,3))
hold off

% CornersL  = [446, 221;
%  466, 225;
%  460, 244;
%  441, 240];
% 
% CornersR  = [516, 225;
%  537, 227;
%  533, 247;
%  513, 244];
% worldPoints = triangulate(CornersL,CornersR,stereoParams)/-1000;
% pointref = [-0.0332984  0.0196809  -0.289665
% -0.0380471  0.0207732  -0.288544
% -0.0369548   0.020522  -0.283672
% -0.0322061  0.0194297  -0.284792];
% 
% 
% 
% coordCorners  =[
%  -0.032947  0.0248755  -0.292459
% -0.0404506  0.0234013  -0.288994
% -0.0370055  0.0152912  -0.280843
% -0.0301032  0.0168377  -0.284378];
% 
% PointRefMarker  =[
% -0.0459644  0.0129777  -0.218604
% -0.0507031  0.0140904  -0.217461
% -0.0495904  0.0138292  -0.212594
% -0.0448517  0.0127164  -0.213737];
% hold on
% % scatter3(pointref(:,1),pointref(:,2),pointref(:,3));
% scatter3(coordCorners(:,1),coordCorners(:,2),coordCorners(:,3));
% scatter3(worldPoints(:,1),worldPoints(:,2),worldPoints(:,3));
% % scatter3(PointRefMarker(:,1),PointRefMarker(:,2),PointRefMarker(:,3));
% % hold off
% % 
% % hold on
% quiver3(-0.0428246455729628,0.0150777530552112,-0.201239326410742,-0.2245,-0.9561,0.1884);
% quiver3(-0.032947,0.0248755,-0.292459,-0.168071,-0.730693,-0.661695);
% quiver3(0,0,0,0,0,-1);
% axis equal
% grid on
% hold off 
% coordCorners  =[
% -0.0321525  0.0200413   -0.28873
% -0.0407938   0.018863  -0.288731
% -0.0388298  0.0110072   -0.28873
% -0.0301581  0.0120168  -0.284732];
% CornersL  =[445, 233;
%  467, 236;
%  462, 256;
%  441, 253];
% 
% CornersR  =[516, 236;
%  538, 239;
%  533, 259;
%  513, 256];
% worldPoints = triangulate(CornersL,CornersR,stereoParams)/-1000;
% 
% norm(coordCorners(1,:) - coordCorners(2,:))
% norm(coordCorners(1,:) - coordCorners(3,:))
% norm(coordCorners(1,:) - coordCorners(4,:))
% 
% norm(worldPoints(1,:) - worldPoints(2,:))
% norm(worldPoints(1,:) - worldPoints(3,:))
% norm(worldPoints(1,:) - worldPoints(4,:))
% hold on
% % scatter3(pointref(:,1),pointref(:,2),pointref(:,3));
% scatter3(coordCorners(:,1),coordCorners(:,2),coordCorners(:,3));
% scatter3(worldPoints(:,1),worldPoints(:,2),worldPoints(:,3));
% axis equal
% grid on
% % scatter3(PointRefMarker(:,1),PointRefMarker(:,2),PointRefMarker(:,3));


% A = rgb2gray(imread('Kidney0900.bmp'));
% A2 = A(100:400,350:830);
% A3 = medfilt2(A,[25 25]);
% J = adapthisteq(A3);
% A_sharp = imsharpen(A3,'Radius',20,'Amount',2);
% A_modified = A_sharp - A3;
% % L = watershed(A_modified);
% BW = imbinarize(A_modified,'global');
% % [centers,radii] = imfindcircles(BW,[50 100],'Sensitivity',0.95);
% imshow(BW);
% % viscircles(centers,radii);
% hblob = vision.BlobAnalysis('AreaOutputPort', true, ... % Set blob analysis handling
%                                 'CentroidOutputPort', true, ... 
%                                 'BoundingBoxOutputPort', true', ...
%                                 'MinimumBlobArea', 5000, ...
%                                 'MaximumBlobArea', 10000, ...
%                                 'MaximumCount', 1);
% 
% [area,centroid, bbox] = hblob(BW);
% bbox = bbox + int32([350,100,0,0]);
% tracked_photo = insertShape(A, 'Rectangle', bbox); 
% 
% s = regionprops(BW,{...
%     'Centroid',...
%     'MajorAxisLength',...
%     'MinorAxisLength',...
%     'Orientation'});
% t = linspace(0,2*pi,50);
% imshow(tracked_photo)
% 
% hold on
% for k = 1:length(s)
%     a = s(k).MajorAxisLength/2;
%     b = s(k).MinorAxisLength/2;
%     Xc = s(k).Centroid(1)+350;
%     Yc = s(k).Centroid(2)+100;
%     phi = deg2rad(-s(k).Orientation);
%     x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
%     y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
%     plot(x,y,'r','Linewidth',5)
% end
% hold off
