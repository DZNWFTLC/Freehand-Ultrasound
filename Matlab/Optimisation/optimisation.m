clear
Mark = rgb2gray(imread('4.png'));
Image = rgb2gray(imread('test3.png'));
Corners_mat = readmatrix('data_corner.csv');
[Num_data,~] = size(Corners_mat);
% figure(1)
% imshow(Mark);
Options = 3;
Scale = 1;
Base_mat = ([0:1:Options-1]-floor(Options/2))*Scale;
Matrix = [];
Similarity = [];
for i = 1:8
    seed = [];
    for j = 1:Options
        seed = [seed; repmat(Base_mat(j),Options^(i-1),1)];
    end
    index = repmat(seed,Options^(8-i),1);
    Matrix = [Matrix index];
end
for num = 1:Num_data
    if num > 99
        filename = ['Kidney0',num2str(num),'.png'];
    elseif num>9
        filename = ['Kidney00',num2str(num),'.png'];
    else
        filename = ['Kidney000',num2str(num),'.png'];
    end
    
    CameraImage= rgb2gray(imread(['CameraFrames\',filename]));
    Corners_mat_temp = repmat(Corners_mat(num,:),length(Matrix),1) + Matrix;

    for i = 1:length(Matrix)
        Coord_Col = Corners_mat_temp(i,1:4);
        Coord_Row = Corners_mat_temp(i,5:end);
        Range = [min(Coord_Row),max(Coord_Row);
            min(Coord_Col),max(Coord_Col)];
        Window = CameraImage(Range(1,1)-10:Range(1,2)+10,Range(2,1)-10:Range(2,2)+10);
        %     figure(2)
        %     imshow(Window);
        
        
        %Window = imbilatfilt(Window);
        Coord = [Coord_Col.'-min(Coord_Col) Coord_Row.'-min(Coord_Row)];
        [Row,Col] = size(Mark);
        Coord_ori = [0 0; Row 0; Row Col;0 Col;];
        tform = estimateGeometricTransform(Coord_ori,Coord,'projective');
        Template = imwarp(imcomplement(Mark),tform, 'SmoothEdges', true);
        Col_valid = any(Template,1);
        Row_valid = any(Template,2);
        Template = imcomplement(Template(Row_valid,Col_valid));
        %     ratio_x =  / rect.width;
        % 	ratio_y = (float)dist_y / rect.height;
        %
        % 	if (abs(ratio_x - 1) <= abs(ratio_y - 1))
        % 		Resize = Size(round(ratio_x * rect.width), round(ratio_x * rect.height));
        % 	else
        % 		Resize = Size(round(ratio_y * rect.width), round(ratio_y * rect.height));
        % 	resize(img_perspective(rect), img_marker, Resize, 0, 0, INTER_AREA);
        
        %imresize
        %     figure(3)
        %     imshow(Template);
        
        c = normxcorr2(Template,Window);
        % figure, surf(c), shading flat
        [ypeak, xpeak] = find(c==max(c(:)));
        
        yoffSet = ypeak-size(Template,1);
        xoffSet = xpeak-size(Template,2);
        
        yoffSet = gather(ypeak-size(Template,1));
        xoffSet = gather(xpeak-size(Template,2));
        Similarity = [Similarity;max(c(:))];
        if max(c(:))>= 0.85
            break;
        end
        % figure
        % imshow(Window);
        % imrect(gca, [xoffSet+1, yoffSet+1, size(Template,2), size(Template,1)]);
    end
end
[Max, index] = max(Similarity)
Corners_mat_temp(index,:);
Corr_ori = Corner_Optimisation(Corners_mat(num,:),CameraImage)
Corr_opt = Corner_Optimisation(Corners_mat_temp(index,:),CameraImage)
