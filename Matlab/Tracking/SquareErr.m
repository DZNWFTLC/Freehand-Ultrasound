function SSE = SquareErr(Valid)

GroundTruth = [0 25 50 25.0147 50.3850;
               25 0 25 5.6311 25.7614;
               50 25 0 26.1485 6.2169;
               25.1047 5.6311 26.1485 0 25.535;
               50.3850 25.7614 6.2169 25.535 0];
Dist = zeros(5);
for i = 1:5
    for j = i+1:5
        dist = norm(Valid(i,:)-Valid(j,:));
        Dist(i,j) = dist;
        Dist(j,i) = dist;
    end
end
SE = (GroundTruth-Dist).^2;
SSE = sum(SE(:));