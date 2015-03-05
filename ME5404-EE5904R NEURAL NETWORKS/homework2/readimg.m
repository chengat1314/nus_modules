function [trainB1,trainB2,testB1,testB2]= readimg()

path_train_glass = 'E:/courses/Course 2014 semester 2/ME5404-EE5904R NEURAL NETWORKS/homework2/Glasswearing/Train/WithGlasses';
path_train_noglass = 'E:/courses/Course 2014 semester 2/ME5404-EE5904R NEURAL NETWORKS/homework2/Glasswearing/Train/WithoutGlasses';
path_test_glass = 'E:/courses/Course 2014 semester 2/ME5404-EE5904R NEURAL NETWORKS/homework2/Glasswearing/Test/WithGlasses';
path_test_noglass = 'E:/courses/Course 2014 semester 2/ME5404-EE5904R NEURAL NETWORKS/homework2/Glasswearing/Test/WithoutGlasses';
dir1 = dir(path_train_glass);
dir2 = dir(path_train_noglass);
dir3 = dir(path_test_glass);
dir4 = dir(path_test_noglass);
for i = 1:length(dir1)-2
    ind = dir1(i+2);
    filename = getfield(ind,'name');
    file =strcat(path_train_glass,'/',filename) ;
    A = double(imread(file));
    trainB1(:,i) = A(:);
%     imshow(A,[]); 
%     pause(2)
end

for i = 1:length(dir2)-2
    ind = dir2(i+2);
    filename = getfield(ind,'name');
    file =strcat(path_train_noglass,'/',filename) ;
    A = double(imread(file));
    trainB2(:,i) = A(:); 
end

for i = 1:length(dir3)-2
    ind = dir3(i+2);
    filename = getfield(ind,'name');
    file =strcat(path_test_glass,'/',filename) ;
    A = double(imread(file));
    testB1(:,i) = A(:); 
end

for i = 1:length(dir4)-2
    ind = dir4(i+2);
    filename = getfield(ind,'name');
    file =strcat(path_test_noglass,'/',filename) ;
    A = double(imread(file));
    testB2(:,i) = A(:); 
end