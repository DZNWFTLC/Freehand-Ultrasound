clear
% RTvector = readmatrix('datasinglenew3.csv');
Cornerdata = readmatrix('datastereonew4.csv');
load('Stereo.mat')
translation = [];
rotation = [];
hold on
for i = 1:length(Cornerdata)
    CornersL = reshape(Cornerdata(i,1:8),[4,2]);
    CornersR = reshape(Cornerdata(i,9:end),[4,2]);
    worldPoints(:,:,i) = triangulate(CornersL,CornersR,stereoParams)/1000;
    centroid(i,:) = mean(worldPoints(:,:,i));
    [p,~,rp] = regress(ones(size(worldPoints(:,:,i),1),1),worldPoints(:,:,i));
    r = vrrotvec(p/norm(p),[0,0,-1]);
        R = vrrotvec2mat(r);

   
%     R_ref = reshape(Rdata(i,:),[3,3]).';
    
    T = [R zeros(3,1);
        0 0 0 1];
    R_dir = T*[0;1;0;0];
    quiver3(centroid(i,1),centroid(i,2),centroid(i,3),R_dir(1),R_dir(2),R_dir(3),0.003);
    
    if (i > 1)
        G = R*R_old.';
        A = (trace(G)-1)/2;
        if A > 1
            A = 1;
        elseif A < -1
            A = -1;
        end
        theta(i-1) = acos(A);
        if theta(i-1)*57.2957795<30
            R_old = R;
        end
        Dist(i-1) = norm(T_old - centroid(i,:));
        if Dist(i-1)<0.06
            T_old = centroid(i,:);
        end
        
    else
        R_old = R;
        T_old = centroid(i,:);
    end
end
plot3(centroid(:,1),centroid(:,2),centroid(:,3));

axis equal
grid on
hold off

% plot3(RTvector(:,1),RTvector(:,2),RTvector(:,3));
% hold on
% axis equal
% grid on
% plot3(centroid(:,1),centroid(:,2),centroid(:,3));
% hold off

max(theta*57.2957795)
plot(theta*57.2957795)
title('Relative Rotation between poses')
xlabel('Frames')
ylabel('Angle/ degree')
Incorrect = length(find(theta*57.2957795>30));
mean(theta)*57.2957795

% max(Dist)
% plot(Dist)
% title('Relative Displacement between Poses')
% xlabel('Frames')
% ylabel('Displacement/ m')
% Incorrect = length(find(Dist>0.006));
% mean(Dist)