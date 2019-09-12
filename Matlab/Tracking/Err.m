function SE = Err(Valid)

GroundTruth = [
    -5.6 0 -2.5;
    -5.6 0 -27.5;
    -5.6 0 -52.5;
    0 -2.7 -52.5;
    0 -2.11 -26.965;
];
SE = 0;
for i = 1:5
SE = SE + norm(Valid(i,:)-GroundTruth(i,:));
end