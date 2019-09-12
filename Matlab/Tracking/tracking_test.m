image = imread('IMG_0205.jpg');
Thresh = 0.05; % Threshold
Thresh2 = 0.1;
hblob = vision.BlobAnalysis('AreaOutputPort', true, ... % Set blob analysis handling
                                'CentroidOutputPort', true, ... 
                                'BoundingBoxOutputPort', true', ...
                                'MinimumBlobArea', 200, ...
                                'MaximumBlobArea', 3000, ...
                                'MaximumCount', 1);

grey = rgb2gray(image);

channelRed = image(:,:,1);
channelGreen = image(:,:,2);
channelBlue = image(:,:,3);

ratio_Y = [0.5 0.5 0];
channelYellow = channelRed*ratio_Y(1) + channelGreen*ratio_Y(2) + channelBlue*ratio_Y(3);

diffFrameGreen = imsubtract(channelGreen, grey); 
figure(1)
imshow(diffFrameGreen);

diffFrameYellow =  imsubtract(channelYellow, grey);
figure(2)
imshow(diffFrameYellow);

diffFrameBlue = imsubtract(channelBlue, grey); 
figure(3)
imshow(diffFrameBlue);
% diffFrameGreen = medfilt2(diffFrameGreen, [3 3]); % Filter out the noise by using median filter

bin_Yellow = imbinarize(diffFrameYellow, Thresh); % Convert the image into binary image with the green objects as white
[~, centroid, bboxYellow] = hblob(bin_Yellow);
tracked_photo = insertShape(image, 'Rectangle', bboxYellow); % Instert the red box
tracked_photo = insertText(tracked_photo,centroid,['X=',num2str(centroid(1)),',Y=',num2str(centroid(2))],'Font','LucidaBrightRegular','BoxColor','w');

bin_Green = imbinarize(diffFrameGreen, Thresh2); % Convert the image into binary image with the green objects as white
[~, centroid, bboxGreen] = hblob(bin_Green);
tracked_photo = insertShape(tracked_photo, 'Rectangle', bboxGreen); % Instert the red box
tracked_photo = insertText(tracked_photo,centroid,['X=',num2str(centroid(1)),',Y=',num2str(centroid(2))],'Font','LucidaBrightRegular','BoxColor','w');

bin_Blue = imbinarize(diffFrameBlue, Thresh); % Convert the image into binary image with the green objects as white
[~, centroid, bboxBlue] = hblob(bin_Blue);
tracked_photo = insertShape(tracked_photo, 'Rectangle', bboxBlue); % Instert the red box
tracked_photo = insertText(tracked_photo,centroid,['X=',num2str(centroid(1)),',Y=',num2str(centroid(2))],'Font','LucidaBrightRegular','BoxColor','w');

figure(4)
imshow(tracked_photo);