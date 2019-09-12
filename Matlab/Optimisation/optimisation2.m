clear
global Window;
Mark = rgb2gray(imread('4.png'));
Image = rgb2gray(imread('test3.png'));
Corners_mat = readmatrix('data.csv');
[Num_data,~] = size(Corners_mat);
% figure(1)
% imshow(Mark);
for num = 1:Num_data
    if num > 99
        filename = ['Kidney0',num2str(num),'.png'];
    elseif num>9
        filename = ['Kidney00',num2str(num),'.png'];
    else
        filename = ['Kidney000',num2str(num),'.png'];
    end
    
    CameraImage= rgb2gray(imread(['CameraFrames\',filename]));
        Coord_Col = Corners_mat(num,1:4);
        Coord_Row = Corners_mat(num,5:end);
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
        fun = @Corner_Optimisation2;
 options = optimoptions('fsolve','Algorithm','levenberg-marquardt', 'Display','off','FiniteDifferenceStepSize',0.001);%,'MaxIterations',1e10,'StepTolerance',1e-10,'MaxFunctionEvaluations',1e10,'FunctionTolerance',1e-10);

        tform_opt = fsolve(fun,tform.T,options)

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
        % figure
        % imshow(Window);
        % imrect(gca, [xoffSet+1, yoffSet+1, size(Template,2), size(Template,1)]);
    
end
