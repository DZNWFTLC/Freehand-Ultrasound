
RTvector = readmatrix('datasingle.csv');
Cornerdata = readmatrix('datastereo.csv');
load('Stereo.mat')
translation = [];
rotation = [];
% axis equal
% grid on
% hold on
% % scatter3(0,0,0);
% for i = 1:length(RTvector)-1
%     Tvec = RTvector(i,1:3).';
%     Rvec = RTvector(i,4:6).';
%     angle = norm(Rvec);
%     Rvec_update = Rvec/angle;
%     R = eye(3)*cos(angle)+(1-cos(angle))*Rvec_update*Rvec_update.'+sin(angle)*[0 -Rvec_update(3) Rvec_update(2);...
%         Rvec_update(3) 0 -Rvec_update(1);...
%         -Rvec_update(2) Rvec_update(1) 0];
%     
% %     R_ref = reshape(Rdata(i,:),[3,3]).';
%     rot_angles = rotm2eul(R);
%     transformation = [Tdata(i,1:4);Tdata(i,5:8);Tdata(i,9:12)];
%     rotation = transformation(:,1:3);
% 
%     T = [rotation Tvec;
%         0 0 0 1];
%     R_dir = T*[0;1;0;0];
%     quiver3(Tvec(1),Tvec(2),Tvec(3),R_dir(1),R_dir(2),R_dir(3),0.003);
%     
%     if (i > 1)
%         G = R*R_old.';
%         A = (trace(G)-1)/2;
%         if A > 1
%             A = 1;
%         elseif A < -1
%             A = -1;
%         end
%         theta(i-1) = acos(A);
%         if theta(i-1)*57.2957795<30
%             R_old = R;
%         end
%         Dist(i-1) = norm(T_old - Tvec);
%         if Dist(i-1)<0.06
%             T_old = Tvec;
%         end
%         
%     else
%         R_old = R;
%         T_old = Tvec;
%     end
% end
% hold off
% max(theta*57.2957795)
% plot(theta*57.2957795)
% title('Relative Rotation between poses')
% xlabel('Frames') 
% ylabel('Angle/ degree') 
% Incorrect = length(find(theta*57.2957795>30));
% max(Dist)
% plot(Dist)
% title('Relative Displacement between Poses')
% xlabel('Frames') 
% ylabel('Displacement/ m') 
% Incorrect = length(find(Dist>0.006));
hold on
for i = 1:length(RTvector)-1
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
    R_dir = T*[0;-1;0;0];
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

plot3(RTvector(:,1),RTvector(:,2),RTvector(:,3));

hold on
axis equal
grid on
plot3(centroid(:,1),centroid(:,2),centroid(:,3));

% hold off
% max(theta*57.2957795)
% plot(theta*57.2957795)
% title('Relative Rotation between poses')
% xlabel('Frames') 
% ylabel('Angle/ degree') 
% Incorrect = length(find(theta*57.2957795>30));
% max(Dist)
% plot(Dist)
% title('Relative Displacement between Poses')
% xlabel('Frames') 
% ylabel('Displacement/ m') 
% Incorrect = length(find(Dist>0.006));