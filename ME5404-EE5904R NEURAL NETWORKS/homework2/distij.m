function dist = distij(data,center)
[m1,n1] = size(data);
[m2,n2] =size(center);
if n1~=n2
    fprintf('dimension different');
    exit;
end
dist= zeros(m1,m2);
for i = 1:m1
    for j = 1:m2
        pointi = data(i,:);
        pointj = center(j,:);
        dis = exp(-sum((pointi-pointj).^2));
        dist(i,j)=dis;
%         dist(j,i)=dis;
    end
end