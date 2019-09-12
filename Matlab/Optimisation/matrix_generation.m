size = 3;
Base_mat = [0:1:size-1]-floor(size/2);
Matrix= [];
for i = 1:8
    seed = [];
    index = [];
    for j = 1:size
    seed = [seed; repmat(Base_mat(j),size^(i-1),1)];
    end
    index = repmat(seed,size^(8-i),1);
    Matrix = [Matrix index];
end