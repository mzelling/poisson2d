function out = eval_sol(X,Y,coef,vert,tria)
%eval_sol Evaluate finite element solution to Laplace's equation
%
%   Z = eval_sol(X,Y,coef,vert,tria) returns a matrix Z whose i,j entry is
%   the finite element solution evaluated at coordinates X(i,j), Y(i,j),
%   where X, Y are matrices of the same size that define a grid (e.g., the 
%   built-in meshgrid function can be used to conveniently define X and Y).
%
%   The argument coef is a row vector giving the solution coefficients 
%   in the finite element function space, as determined by laplace2d. Thus,
%   coef has n_nodes rows, where n_nodes is the number of vertices in the 
%   discretization. The other arguments specify the vertices and triangles 
%   of the discretization; vert is an n_nodes-by-2 array that lists the 
%   coordinates for each point in the discretization; tria is an 
%   n_triangles-by-3 array whose i-th row specifies the (global) indices of
%   the three vertices in the i-th triangle of the discretization.

out = -13*ones(size(X));

for j = 1:size(X,2)
    % Construct vector of points to evaluate
    pts = [X(:,j), Y(:,j)];
    
    % Identify the triangles containing the points
    [vert2triangle,triangle_idx] = findtria(vert,tria,pts);
    % For more detail on the preceding two lines, see comment at the bottom
    % of this file.

    for i = 1:size(pts,1)
        % Check if point lies in no triangle
        if vert2triangle(i,1) == 0
            out(i,j) = 0; % apply Dirichlet boundary condition
        else
            % Transform mesh coordinates of point i to the reference triangle
            local2global = tria(triangle_idx(vert2triangle(i,1)),:);
            this_pt = pts(i,:);
            J = [vert(local2global(2),:)-vert(local2global(1),:);
                 vert(local2global(3),:)-vert(local2global(1),:)]';
            ref_coords = J\((this_pt - vert(local2global(1),:))');

            % Evaluate reference basis functions at this point
            ref_basis_1 = 1 - ref_coords(1) - ref_coords(2);
            ref_basis_2 = ref_coords(1);
            ref_basis_3 = ref_coords(2);

            % Compute solution value
            out(i,j) = coef(local2global(1))*ref_basis_1 + ...
                       coef(local2global(2))*ref_basis_2 + ...
                       coef(local2global(3))*ref_basis_3;
        end
    end
end

% Comment on triangle-finding:
%   Since the documentation for the findtria function in mesh2d's aabb-tree 
% package is difficult to understand at first reading, we explain the code
% here. 
%   The vector triangle_idx contains a vector that list the triangles
% each queried point may or may not be contained in. In case a point lies
% at the intersection of several triangles, triangle_idx contains multiple
% successive entries corresponding to this point. 
%   The array 
% vert_to_triangle associates each vertex to all the triangles it
% intersects; it has nrow(pts) rows and 2 columns. 
%   For evaluating a finite element solution, we only need one triangle 
% intersecting each queried point, so we select triangle indices with the 
% first column of vert_to_triangle.

