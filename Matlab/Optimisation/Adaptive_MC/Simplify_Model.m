function [model_new] = Simplify_Model(model,Cube_size)
% SIMPLIFY_MODEL - Merge large CT scan model into smaller one
%
% Inputs:
%    model - 3D CT images matrix
%    Cube_size - minimum size of a cube
%
% Outputs:
%    model_new - updated 3D CT images matrix
%
%
% Other m-files required: run.m
% Subfunctions: Intersection_Possion, Face_Vertex
% MAT-files required: none
%
% Author: Chongyun WANG
% email: chongyun.wang.18@ucl.ac.uk
% Last revision: May-2019
%------------- BEGIN CODE --------------
[I,J,K] = size(model);

if logical(mod(I,Cube_size))
    I_add = Cube_size-mod(I,Cube_size);
else
    I_add = 0;
    
end

if logical(mod(J,Cube_size))
    J_add = Cube_size-mod(J,Cube_size);
else
    J_add = 0;
end

if logical(mod(K,Cube_size))
    K_add = Cube_size-mod(K,Cube_size);
else
    K_add = 0;
end

model = padarray(model,[I_add,J_add,K_add],'symmetric','post');
[I,J,K] = size(model);

SUM=zeros(I+1-Cube_size,J+1-Cube_size,K+1-Cube_size);
for i = 1:Cube_size
    for j = 1:Cube_size
        for k = 1:Cube_size
            SUM = SUM + model(i:end-Cube_size+i,j:end-Cube_size+j,k:end-Cube_size+k);
        end
    end
end
[I,J,K] = size(model);

x = [];
num = 1;
for i = 1:Cube_size:I
    x(num,:,:) = SUM(i,:,:);
    num = num + 1;
end

y = [];
num = 1;
for j = 1:Cube_size:J
    y(:,num,:) = x(:,j,:);
    num = num +1;
end

model_new = [];
num = 1;
for k = 1:Cube_size:K
    model_new(:,:,num) = y(:,:,k);
    num = num + 1;
end

model_new = model_new/Cube_size^3;