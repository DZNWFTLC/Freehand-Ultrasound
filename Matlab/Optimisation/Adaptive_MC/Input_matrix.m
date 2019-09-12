function [Cube_map] = Input_matrix(model,isovalue,isPoisson)
% INPUT_MATRIX - Perform marching cubes algorithm using output
% from possion algorithm
%
% Inputs:
%    model - 3 dimensional CT images matrix
%    isovalue - threshold for identifying whether vertices are inside the
%    surface or not
%    isPoisson - logical flag
% Outputs:
%    Cube_map - 3 dimensional binary matrix
%
%
% Other m-files required: run.m
% Subfunctions: none
% MAT-files required: none
%
% Author: Chongyun WANG
% email: chongyun.wang.18@ucl.ac.uk
% Last revision: May-2019
%------------- BEGIN CODE --------------


Cube_num = size(model) - 1; % number of cubes along each direction of image

Cube_map = zeros(Cube_num(1),Cube_num(2),Cube_num(3),'uint16');

% Cube vertices
%           5           6   Z
%       8           7
%       
%       
%
%           1           2   X
%       4           3
%       Y

vertex_index = {1:Cube_num(1), 1:Cube_num(2), 1:Cube_num(3); ... % create index matrix for vertices of each cube
    2:Cube_num(1)+1, 1:Cube_num(2), 1:Cube_num(3); ...
    2:Cube_num(1)+1, 2:Cube_num(2)+1, 1:Cube_num(3); ...
    1:Cube_num(1), 2:Cube_num(2)+1, 1:Cube_num(3); ...
    1:Cube_num(1), 1:Cube_num(2), 2:Cube_num(3)+1; ...
    2:Cube_num(1)+1, 1:Cube_num(2), 2:Cube_num(3)+1; ...
    2:Cube_num(1)+1, 2:Cube_num(2)+1, 2:Cube_num(3)+1; ...
    1:Cube_num(1), 2:Cube_num(2)+1, 2:Cube_num(3)+1 };


for i= 1:8                           
    if isPoisson
        inside = model(vertex_index{i, :}) == 1;
    else
        inside = model(vertex_index{i, :}) > isovalue;  % check if vertices are inside the surface
    end
    
    Cube_map(inside) = bitset(Cube_map(inside), i);     % generate Cube matrix ( 8bit for each cube)
end
