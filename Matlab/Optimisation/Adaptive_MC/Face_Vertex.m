function [Face, Vertices] = Face_Vertex(P,Point_index,List_triangle)
Face = [];
Vertices = [];

for i=1:3:13
    Triangle_valid_index = find(List_triangle(:, i)>-1);
    % find valid triangle surface
    
    Vertices = List_triangle(Triangle_valid_index, i:i+2)+1;
    % every three element in List_triangle corresponding to vertices of a
    % triangle 
    
    if ( ~ isempty(Triangle_valid_index) )
        Vertex1_index = sub2ind(size(Point_index), Triangle_valid_index, Vertices(:,1));
        Vertex2_index = sub2ind(size(Point_index), Triangle_valid_index, Vertices(:,2));
        Vertex3_index = sub2ind(size(Point_index), Triangle_valid_index, Vertices(:,3));
        % get index of vertices 
        
        Face = [Face; [Point_index(Vertex1_index) Point_index(Vertex2_index) Point_index(Vertex3_index)]];
    end
    
end

for i = 1:12
    Vertex_index = any(Point_index(:, i),2);
    if (~ isempty(Vertex_index))
        Vertices(Point_index(Vertex_index, i), 1:3) = P(Vertex_index, :, i);
    end
end