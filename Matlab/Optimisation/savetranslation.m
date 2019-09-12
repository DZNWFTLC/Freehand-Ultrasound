%Kinematic model of da vinci

%opens file and formats into matrix with 12 columns
fileID = fopen('test.txt');
sizeA = [36 Inf];
formatSpec = '%f';
A = fscanf(fileID,formatSpec,sizeA);
A = A';

%deletes extra columns
V = 13:1:36;
A(:,V) = [];

%position vector
v = [0 0 0 1]';

%file to write to
fileID = fopen('transform.txt','w');

len = length(A(:,1));
for i=1:len
    
    %create a 4x4 matrix with last row 0001
    M = A(i,1:end);
    M = reshape(M,4,3);
    M = M';
    insert = [0 0 0 1];
    
    %T is the transformation matrix
    T = [M; insert];
    
    position = T*v;
    x = position(1);
    y = position(2);
    z = position(3);
     
    %test for transformation
    rot = -4.35;
    translation = [0.01231 0 0.00762]';
 
    %makes a transformation matrix
    Rotation = roty(rot);
    T_g2p = [Rotation translation; v'];
    
    %both transformations combined
    Tcombined = T*R;
    ultra = Tcombined*v;
    
    line = [Tcombined(1,:) Tcombined(2,:) Tcombined(3,:) 0 0 0 1];
    
    fprintf(fileID,'%6.5f ',line);
    fprintf(fileID, '\n');
end

fclose(fileID);