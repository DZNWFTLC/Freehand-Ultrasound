Valid =[
    5.5384    0.7026  119.8596;
    21.6027   12.8981  108.9038;
    36.8732   23.1884   95.1366;
    39.8395   16.3060   92.0798;
    25.4536    6.7102  111.5937;
];

GroundTruth = [
    -5.6 0 -2.5;
    -5.6 0 -27.5;
    -5.6 0 -52.5;
    0 -2.7 -52.5;
    0 -2.11 -26.965;
];

PC_V = pointCloud(Valid);
PC_GT = pointCloud(GroundTruth);

% tform = pcregistericp(PC_V,PC_GT,'MaxIterations',200,'Tolerance',[0.001 0.005]);
tform = pcregistercpd(PC_V,PC_GT,'Transform','Rigid','MaxIterations',200,'Tolerance',1e-20);

invtform = invert(tform);
ptCloudOut = pctransform(PC_GT,invtform);
pcshow(ptCloudOut,'MarkerSize', 100)
    grid on
    axis equal
    hold on 
    pcshow(PC_V,'MarkerSize', 500)
hold off

% plot3(Valid(:,1),Valid(:,2),Valid(:,3));
% plot3(GroundTruth(:,1),GroundTruth(:,2),GroundTruth(:,3));
%     axis equal
%     grid on