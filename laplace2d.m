function u_h = laplace2d(varargin)
% laplace2d Solve Laplace's equation on 2D polygonal domain
%
% u_h = laplace2d(polygon_vert,f,h,polygon_edges) returns a function handle
% that evaluates the finite element solution for the PDE -Laplace(u) = f 
% subject to Dirichlet boundary conditions (i.e. u_h evaluates to zero on 
% the polygon's boundary).
%
% The argument polygon_vert supplies the polygonal vertices as a k-by-2 
% array, where k is the number of vertices and  polygon_vert(i,:) gives 
% the 2D coordinates of the i-th vertex. 
%
% The argument f supplies the forcing function (right hand side of PDE) as 
% a function handle of the form @(x,y) (...), where x,y are the xy
% coordinates for a single point in 2D.
%
% The argument h controls the fineness of the discretized grid. More
% precisely, it is an upper bound on the length of triangular edges in the 
% discretization.
%
% The argument polygon_edges is an array with two columns whose i-th row 
% specifies the indices in [1,2,3,...,k] of the nodes connected by the i-th
% edge. If this argument is omitted, it is assumed that the nodes are 
% connected in ascending order.

% Set arguments
polygon_vert = varargin{1};
f = varargin{2};
h = varargin{3};

m = size(polygon_vert,1);
if (nargin == 4)
    polygon_edges = varargin{4};
else
    polygon_edges = repmat([1,2],m,1) + (0:(m-1))';
    polygon_edges(m,2) = 1;
end

% Construct domain triangularization
[vert,~,triangles,~] = refine2(polygon_vert,polygon_edges,[],[],h);
n_triangles = size(triangles,1);
n = size(vert,1);

% Identify boundary nodes by creating an enclosing polygon. If a point is
% both in the original polygon and the enclosing polygon, it is a boundary
% point
M = 2*max(polygon_vert,[],'all');
vert_out = [polygon_vert;-M,-M; M,-M; M,M; -M,M];
edges_out = [polygon_edges; m+1,m+2; m+2,m+3; m+3,m+4; m+4,m+1];
bdry_idx = inpoly2(vert,vert_out,edges_out);

% Pre-allocate sparse stiffness matrix
A = sparse(n,n);
b = zeros(n,1);

for K = 1:n_triangles
    % specify vertices of this triangle
    loc2glob = triangles(K,:); % local-to-global index mapping
    vert_loc = [vert(loc2glob(1),:); 
                vert(loc2glob(2),:);
                vert(loc2glob(3),:)];
    % compute local tensor, including local forcing
    [dA, db] = local_tensor(vert_loc,f);
    
    % update stiffness matrix
    A(loc2glob, loc2glob) = A(loc2glob, loc2glob) + dA; 
    b(loc2glob) = b(loc2glob) + db;
end

% Impose boundary conditions and solve Ax = b
coef = zeros(n,1);
coef_inside = A(~bdry_idx,~bdry_idx)\b(~bdry_idx);
coef(~bdry_idx) = coef_inside;

% Compute function handle to evaluate PDE solution
u_h = @(X,Y) eval_sol(X,Y,coef,vert,triangles);

end
