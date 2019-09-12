clear
ImageL = imread('53779L.png');
ImageR = imread('53779R.png');
Image = {ImageL(46:end,1:end-5,:),ImageR(46:end,1:end-5,:)};
load('Stereo.mat');
Thresh = 0.05; % Threshold
Thresh2 = 0.1;

GroundTruth = [0 25 50 25.0147 50.3850;
               25 0 25 5.6311 25.7614;
               50 25 0 26.1485 6.2169;
               25.1047 5.6311 26.1485 0 25.535;
               50.3850 25.7614 6.2169 25.535 0];
           
           
% hblob = vision.BlobAnalysis('AreaOutputPort', true, ... % Set blob analysis handling
%     'CentroidOutputPort', true, ...
%     'BoundingBoxOutputPort', true', ...
%     'MinimumBlobArea', 200, ...
%     'MaximumBlobArea', 500, ...
%     'MaximumCount', 1);


for i = 1:2
    image = Image{i};
    Centroid = [];
    grey = rgb2gray(image);
    
    channelRed = image(:,:,1);
    channelGreen = image(:,:,2);
    channelBlue = image(:,:,3);
    
    ratio_Y = [0.5 0.5 0];
    channelYellow = channelRed*ratio_Y(1) + channelGreen*ratio_Y(2) + channelBlue*ratio_Y(3);
    
    diffFrameGreen = imsubtract(channelGreen, grey);
    % figure(1)
    % imshow(diffFrameGreen);
    
    diffFrameYellow =  imsubtract(channelYellow, grey);
    % figure(2)
    % imshow(diffFrameYellow);
    
    diffFrameBlue = imsubtract(channelBlue, grey);
    % figure(3)
    % imshow(diffFrameBlue);
    % diffFrameGreen = medfilt2(diffFrameGreen, [3 3]); % Filter out the noise by using median filter
    
    bin_Yellow = imbinarize(diffFrameYellow, Thresh); % Convert the image into binary image with the green objects as white
    %      [~, Centroid(1,:,i), bboxYellow] = hblob(bin_Yellow);
    %     tracked_photo = insertShape(image, 'Rectangle', bboxYellow); % Instert the red box
    stats = regionprops(bin_Yellow,'centroid');
    centroid = cell2mat(struct2cell(stats));
    Centroid = [Centroid;reshape(centroid, 2, length(centroid)/2).'];
    
    bin_Green = imbinarize(diffFrameGreen, Thresh2); % Convert the image into binary image with the green objects as white
    %     [~, centroid(2,:,i), bboxGreen] = hblob(bin_Green);
    %     tracked_photo = insertShape(tracked_photo, 'Rectangle', bboxGreen); % Instert the red box
    stats = regionprops(bin_Green,'centroid');
    centroid = cell2mat(struct2cell(stats));
    Centroid = [Centroid;reshape(centroid, 2, length(centroid)/2).'];
    bin_Blue = imbinarize(diffFrameBlue, Thresh2); % Convert the image into binary image with the green objects as white
    %     [~, centroid(3,:,i), bboxBlue] = hblob(bin_Blue);
    %     tracked_photo = insertShape(tracked_photo, 'Rectangle', bboxBlue); % Instert the red box
    stats = regionprops(bin_Blue,'centroid');
    centroid = cell2mat(struct2cell(stats));
    Centroid = [Centroid;reshape(centroid, 2, length(centroid)/2).'];
    
    figure(i)
    imshow(image);
    hold on
    plot(Centroid(:,1),Centroid(:,2),'.');
    hold off
    Centroid_tot{i} = Centroid;
end
if length(Centroid_tot{1})== length(Centroid_tot{2})
    CentroidL = Centroid_tot{1};
    CentroidR = Centroid_tot{2};
    CentroidL(:,1) = CentroidL(:,1) + 45;
    CentroidR(:,1) = CentroidR(:,1) + 45;
    worldPoints = triangulate(CentroidL,CentroidR,stereoParams);
    [L,~] = size(worldPoints);
Valid = [];
for i = 1:L
    if worldPoints(i,3) > 80
        Valid = [Valid; worldPoints(i,:)];
    end
end
    figure(3)
    scatter3(Valid(:,1),-Valid(:,2),Valid(:,3));
    grid on
    axis square
    zlim([0 inf])
    xlim([-50 50])
    ylim([-50 50])
    
end


