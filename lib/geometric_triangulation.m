function [target_nav] = geometric_triangulation(RDs,OCs)
% Function for the geometric solving of triangulation, which gets in input:
% - the ordered list of ray directions in the reference frame (matrix)
% - the ordered list of optical centers in the reference frame (matrix)

A = zeros(3,3);
b = zeros(3,1);
n_pts = size(RDs,1);

for k = 1:n_pts
    % Build the orthogonal projection matrix. Formula: (I - u*u')
    proj_matrix = eye(3) - (RDs(:,k) * RDs(:,k)');
    
    % Sum to obtain the A, b of the linear system to be solved
    A = A + proj_matrix;
    b = b + (proj_matrix * OCs(:,k));
end

% Solve for the Target Position in Navigation Frame (A*x = b)
target_nav = A \ b;


end
