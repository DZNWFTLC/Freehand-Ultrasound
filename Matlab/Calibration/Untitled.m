points = [0.0030838, 0.0211164, -0.365029;
0.018812, 0.0211521, -0.353187;
0.0184867, 0.0378261, -0.365746;
0.00308376, 0.0370072, -0.36503;
0, 0, 0];

hold on
scatter3(points(1,1),points(1,2),points(1,3),'r');
scatter3(points(2,1),points(2,2),points(2,3),'g');
scatter3(points(3,1),points(3,2),points(3,3),'b');
scatter3(points(4,1),points(4,2),points(4,3),'y');
scatter3(points(5,1),points(5,2),points(5,3),'c');
axis equal
grid on
