function [P, Point_index] = Intersection_Poisson(Cube_size,Cube_valid_index,Valid_edge_binary,Vertex_offset,Edge_index)
% INTERSECTION_POSSION - Calculate coordinates of all intersected points
%
% Inputs:
%    Cube_size - size of input matrix
%    Cube_valid_index - Index of nonempty cubes
%    Valid_edge_binary - Logical Index of nonempty cubes
%    Vertex_offset - offset of 8 vertices within a cubd
%    Edge_index - index of 12 edges within a cube
%
% Outputs:
%    P - n*3*12 martrix containing coordinates of intersected points
%    Point_index - n*12 matrix containing index of intersected points
%
%
% Other m-files required: run.m
% Subfunctions: None
% MAT-files required: none
%
% Author: Chongyun WANG
% email: chongyun.wang.18@ucl.ac.uk
% Last revision: May-2019
%------------- BEGIN CODE --------------
P = zeros(length(Cube_valid_index),3,12);
Point_index = zeros(length(Cube_valid_index),12);
num = 0;
for i = 1:12
    Edge_valid_index = Cube_valid_index(Valid_edge_binary(:,i));
    % find valid edges
    [I, J, K] = ind2sub(Cube_size,Edge_valid_index);
    Vertices1 = [I J K] + Vertex_offset(Edge_index(i,1),:);
    Vertices2 = [I J K] + Vertex_offset(Edge_index(i,2),:);
    % coordinates of vertices corresponding to valid edges
    P(Valid_edge_binary(:,i),:,i) = (Vertices1 + Vertices2)/2;
    % coordinates intersection points (for Poisson, take middle point)
    current_index = 1:length(Edge_valid_index);
    
    Point_index(Valid_edge_binary(:,i),i) =  current_index + num;
    % store index of edges
    num = num + length(Edge_valid_index);
end