function [dA,db] = local_tensor(vert_loc, f)
%local_tensor Compute local tensor and forcing for a triangle
%
%   A = local_tensor(vert_loc) returns a 3-by-3 matrix A and 3-by-1 vector 
%   b. The entry A(i,j) is this triangle's contribution to the energy inner 
%   product of the nodal basis functions for vertices i and j. The entry 
%   b(i) is this triangle's contribution to the inner product of the 
%   forcing function f (right hand side of PDE) with the nodal basis 
%   function for the i-th vertex. 
%
%   The vertices vert_loc are specified as a 3-by-2 array whose i-th row 
%   are the coordinates of the i-th vertex in the triangle. 
%   The forcing function f is a function handle of the form @(x,y) (...),
%   where x,y are the xy coordinates of a single point in 2D.

% To evaluate the integral for general vertices, we transform our
% triangle to the reference triangle (the standard two-simplex).

% Jacobian of the transformation to reference triangle.
J = [vert_loc(2,:)-vert_loc(1,:);
    vert_loc(3,:)-vert_loc(1,:)]';
J_scale = abs(det(J)); % scaling factor for transformed integrals
% Gradients of reference nodal basis functions.
ref_gradients = [-1,-1;1,0;0,1]';

dA = zeros(3,3);
db = zeros(3,1);
% Compute energy inner product of nodal basis functions
% for each pair of mesh vertices i,j over this triangle.
% By symmetry, fill only upper triangular part of local tensor.
% Also compute inner product of f with the nodal basis function for
% each vertex.

for i = [1,2,3]
    for j = i:3
        
        % Transformed integrand over reference triangle
        integrand = J_scale*(J'\ref_gradients(:,i))'*(J'\ref_gradients(:,j));

        % Note: integrand is constant over the triangle because the
        % Jacobian and reference gradients are
        
        % Quadrature is easy because integrand is constant
        dA(i,j) = integrand * 0.5; % multiply by area of ref triangle
    end
    
    % To integrate f over triangle, use its value at the midpoint
    triangle_midpt = mean([vert_loc(1,:); vert_loc(2,:); vert_loc(3,:)]);
    
    % Alternatively, average value of f over the vertices
    %integrand = J_scale*f(vert_loc(i,1),vert_loc(i,2))*(1/3);
    
    % Basis functions take value (1/3) at midpoint
    integrand = J_scale*f(triangle_midpt(1),triangle_midpt(2))*(1/3);
    db(i) = integrand * 0.5; % multiply by area of ref triangle
end

dA = dA + triu(dA,+1)';
end

