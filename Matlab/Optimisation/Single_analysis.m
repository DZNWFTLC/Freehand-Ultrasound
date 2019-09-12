clear
RTvector = readmatrix('data.csv');
% Rdata = readmatrix('data2.csv');

plot3(RTvector(:,1),RTvector(:,2),RTvector(:,3));
axis equal
grid on
hold on
% scatter3(0,0,0);
Rotation(:,:,1) = eye(3);
Rotation(:,:,2) = rotz(90);
Rotation(:,:,3) = rotz(180);
Rotation(:,:,4) = rotz(270);
for i = 1:133
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
        [~,I] = min(theta_temp)
        R = Rotation(:,:,I)*R;
        theta(i-1) = acos(I);
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
    
    T = [R zeros(3,1);
        0 0 0 1];
    R_dir = T*[0;-1;0;0];
    quiver3(Tvec(1),Tvec(2),Tvec(3),R_dir(1),R_dir(2),R_dir(3),0.003);
end
hold off
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
