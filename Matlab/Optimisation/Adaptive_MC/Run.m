clear
dimension = 512; % dimension of CT image
num = 192; % number of CT image slices
scale = 1; % scale factor along Z axis
isovalue = 1200; % threshold
Cube_size = 1; % target cube size (for simplify only)

isPoisson = false; % if input comes from Poisson reconstruction
isSimplify = false; % if simplify CT model
ifMean = false; % if apply mean filter

if ~isPoisson
    for i = 1:num
        CT_slide = fopen(strcat('AXIAL\',num2str(i)));
        model(:,:,i) = fread(CT_slide,[dimension,dimension],'uint16','ieee-be');
        fclose(CT_slide);
        % read each CT image sequently and store in a matrix
    end
    if ifMean
        model = medfilt3(model,[3 3 3]);
        % apply 3D mean filter
    end
    
    if isSimplify
        [model] = Simplify_Model(model,Cube_size);
        % apply larger cube size to model
    end
end
Cube_map = Input_matrix(model,isovalue,isPoisson);
% convert CT model to cube matrix
[Face,Vertices] = Marching_Cubes_Poisson(Cube_map,scale);

trimesh(Face,Vertices(:,1),Vertices(:,2),Vertices(:,3))
