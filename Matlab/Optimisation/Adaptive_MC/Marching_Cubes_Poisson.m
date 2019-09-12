function [Face,Vertices] = Marching_Cubes_Poisson(Cube_map,scale)
% MARCHING_CUBES_POSSION - Perform marching cubes algorithm using output
% from possion algorithm
%
% Inputs:
%    Cube_map - 3 dimensional binary matrix
%    scale - scale factor of z axis
%
% Outputs:
%    Face - 3*n matrix containing indices of vertices 
%    Vertices - 3*m matrix containing coordinate of each vertices
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

% Cube vertices
%         L8  5     L5      6   Z
%       8      L7     7   L6
%
%             L9            L10
%       L12           L11
%         L4  1      L1     2   X
%       4      L3     3   L2
%       Y

Edge_index = [1 2;
    2 3;
    3 4;
    4 1;
    5 6;
    6 7;
    7 8;
    8 1;
    1 5;
    2 6;
    3 7;
    4 8];
% edge index defined above

Vertex_offset = [0 0 0;
    1 0 0;
    1 1 0;
    0 1 0;
    0 0 1;
    1 0 1;
    1 1 1;
    0 1 1];
% offset of vertices defined above

Tables;
% get edge table and tri table

Cube_edge = edgeTable(Cube_map+1);  % Check edge table

Cube_size = size(Cube_map);

Cube_valid_index =  find(Cube_edge); % Adaptive MC - find nonzero cubes

% Cube_valid_index =  1:(Cube_size(1)*Cube_size(2)*Cube_size(3)); % Original MC
% Cube_valid_index = Cube_valid_index'; % Original MC

Valid_cube = Cube_map(Cube_valid_index); % find values of nonempty cubes 

Valid_edge = Cube_edge(Cube_valid_index); % find edges corresponding to nonempty cubes

Valid_edge_binary= logical(decimalToBinaryVector(Valid_edge,12,'LSBFirst')); % convert edge values back to binary

[P, Point_index] = Intersection_Poisson(Cube_size,Cube_valid_index,Valid_edge_binary,Vertex_offset,Edge_index);
% find intersected points

List_triangle = triTable(Valid_cube+1,1:15); % for each valid cube, find corresponding triagles
[Face, Vertices] = Face_Vertex(P,Point_index,List_triangle);


[Vertices, Sort_index] = sortrows(Vertices);
[Vertices,~,ic]  = unique(Vertices,'rows');
Vertices = Vertices.*[1; 1; -scale].';
Sort_index(Sort_index) = ic;
Face = unique(Sort_index(Face),'rows');
% find unique vertices
%------------- END CODE --------------