[X,Y] = meshgrid(-1.5:0.05:1.5);
Z = zeros(size(X));
[m,n] = size(Z);
for i=1:m
    for j=1:n
        Z(i,j) = u_h(X(i,j),Y(i,j));
    end
end

surf(X,Y,Z);

