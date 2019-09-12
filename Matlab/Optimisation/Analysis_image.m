

A = rgb2gray(imread('tumour2.png'));
A2 = A(100:400,350:830);
A3 = medfilt2(A2,[25 25]);
J = adapthisteq(A3);
A_sharp = imsharpen(A3,'Radius',20,'Amount',2);
A_modified = A_sharp - A3;
% L = watershed(A_modified);
BW = imbinarize(A_modified,'global');
% [centers,radii] = imfindcircles(BW,[50 100],'Sensitivity',0.95);
imshow(BW);
% viscircles(centers,radii);
hblob = vision.BlobAnalysis('AreaOutputPort', true, ... % Set blob analysis handling
                                'CentroidOutputPort', true, ... 
                                'BoundingBoxOutputPort', true', ...
                                'MinimumBlobArea', 5000, ...
                                'MaximumBlobArea', 10000, ...
                                'MaximumCount', 1);

[area,centroid, bbox] = hblob(BW);
bbox = bbox + int32([350,100,0,0]);
tracked_photo = insertShape(A, 'Rectangle', bbox); 

s = regionprops(BW,{...
    'Centroid',...
    'MajorAxisLength',...
    'MinorAxisLength',...
    'Orientation'});
t = linspace(0,2*pi,50);
imshow(tracked_photo)
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
