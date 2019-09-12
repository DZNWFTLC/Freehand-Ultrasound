rot = -4.35;
translation = [0.01231 0 0.00762]';
load('kinematics.mat')
%makes a transformation matrix
Rotation = rotx(180)*roty(180)*roty(rot);
T_g2p = [Rotation translation; 0 0 0 1];
for i = 1:150
    T = [reshape(A(i,:),[4,3]).';0 0 0 1];
    T_probe = T*T_g2p;
    
    Traj_k(i,:) = T_probe(1:3,4).'+[0.0056, 0.003, -0.0025];
end
plot3(Traj_k(:,2),Traj_k(:,1),Traj_k(:,3));
axis equal
grid on