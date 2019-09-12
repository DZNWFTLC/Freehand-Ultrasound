clear
RTvector = readmatrix('datanew3.csv');
RTvector2 = readmatrix('datanew4.csv');

axis equal
grid on
hold on
scatter3(0,0,0);
Rotation(:,:,1) = eye(3);
Rotation(:,:,2) = rotz(90);
Rotation(:,:,3) = rotz(180);
Rotation(:,:,4) = rotz(270);
for i = 1:length(RTvector)-1
    Tvec = RTvector(i,1:3).';
    Rvec = RTvector(i,4:6).';
    angle = norm(Rvec);
    Rvec_update = Rvec/angle;
    R = eye(3)*cos(angle)+(1-cos(angle))*Rvec_update*Rvec_update.'+sin(angle)*[0 -Rvec_update(3) Rvec_update(2);...
        Rvec_update(3) 0 -Rvec_update(1);...
        -Rvec_update(2) Rvec_update(1) 0];
    
    %     R_ref = reshape(Rdata(i,:),[3,3]).';
    
    
    
    if (i > 1)
        for j = 1:4
            R_temp = R*Rotation(:,:,j);
            G = R_temp*R_old.';
            A = (trace(G)-1)/2;
            if A > 1
                A = 1;
            elseif A < -1
                A = -1;
            end
            theta_temp(j) = acos(A);
        end
        [~,I] = min(theta_temp);
        R = Rotation(:,:,I)*R;
        theta(i-1) = theta_temp(I);
        if theta(i-1)*57.2957795<30
            R_old = R;
        end
        Dist(i-1) = norm(T_old - Tvec);
        if Dist(i-1)<0.06
            T_old = Tvec;
        end
    else
        R_old = R;
        T_old = Tvec;
    end
    
    T = [R Tvec;
        0 0 0 1];
    T_trans = eye(4);
    T_trans(1,4) = -0.025;
    T = T*T_trans;
    R_dir = T*[0;-1;0;0];
    Tvec_new(i,:) = [T(1,4),T(2,4),T(3,4)];

    quiver3(T(1,4),T(2,4),T(3,4),R_dir(1),R_dir(2),R_dir(3),0.003);
end
plot3(Tvec_new(:,1),Tvec_new(:,2),Tvec_new(:,3));
hold off
max(theta*57.2957795)
plot(theta*57.2957795)
title('Relative Rotation between poses')
xlabel('Frames')
ylabel('Angle/ degree')
Incorrect = length(find(theta*57.2957795>30));
mean(theta)*57.2957795

max(Dist)
plot(Dist)
title('Relative Displacement between Poses')
xlabel('Frames')
ylabel('Displacement/ m')
Incorrect = length(find(Dist>0.006));
mean(Dist)

plot3(RTvector2(:,1),RTvector2(:,2),RTvector2(:,3));
hold on
axis equal
grid on
scatter3(0,0,0);
Rotation2(:,:,1) = eye(3);
Rotation2(:,:,2) = rotz(90);
Rotation2(:,:,3) = rotz(180);
Rotation2(:,:,4) = rotz(270);
for i = 1:length(RTvector2)-1
    Tvec2 = RTvector2(i,1:3).';
    Rvec2 = RTvector2(i,4:6).';
    angle2 = norm(Rvec2);
    Rvec_update2 = Rvec2/angle2;
    R2 = eye(3)*cos(angle2)+(1-cos(angle2))*Rvec_update2*Rvec_update2.'+sin(angle2)*[0 -Rvec_update2(3) Rvec_update2(2);...
        Rvec_update2(3) 0 -Rvec_update2(1);...
        -Rvec_update2(2) Rvec_update2(1) 0];
    
    %     R_ref = reshape(Rdata(i,:),[3,3]).';
    
    
    
    if (i > 1)
        for j = 1:4
            R_temp2 = R2*Rotation2(:,:,j);
            G2 = R_temp2*R_old2.';
            A2 = (trace(G2)-1)/2;
            if A2 > 1
                A2 = 1;
            elseif A2 < -1
                A2 = -1;
            end
            theta_temp2(j) = acos(A2);
        end
        [~,I] = min(theta_temp2);
        R2 = Rotation2(:,:,I)*R2;
        theta2(i-1) = theta_temp2(I);
        if theta2(i-1)*57.2957795<30
            R_old2 = R2;
        end
        Dist2(i-1) = norm(T_old2 - Tvec2);
        if Dist2(i-1)<0.06
            T_old2 = Tvec2;
        end
    else
        R_old2 = R2;
        T_old2 = Tvec2;
    end
    
    T2 = [R2 zeros(3,1);
        0 0 0 1];
    R_dir2 = T2*[0;-1;0;0];
    quiver3(Tvec2(1),Tvec2(2),Tvec2(3),R_dir2(1),R_dir2(2),R_dir2(3),0.003);
end
hold off

max(theta2*57.2957795)
plot(theta2*57.2957795)
title('Relative Rotation between poses')
xlabel('Frames')
ylabel('Angle/ degree')
Incorrect2 = length(find(theta2*57.2957795>30));
mean(theta2)*57.2957795

% max(Dist2)
% plot(Dist2)
% title('Relative Displacement between Poses')
% xlabel('Frames')
% ylabel('Displacement/ m')
% Incorrect2 = length(find(Dist2>0.006));
% mean(Dist2)
