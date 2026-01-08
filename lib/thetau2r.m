function R = thetau2r(thetau)

%fprintf(1, 'r2thetau');

theta = norm(thetau);

if(theta < 1e-8)
	R = eye(3);

	return;
end
 
n = thetau/theta;

c = cos(theta);
s = sin(theta);

nx = n(1);
ny = n(2);
nz = n(3);

R = [	c+(1-c)*nx^2		-nz*s+(1-c)*nx*ny	ny*s+(1-c)*nx*nz; 
		nz*s+(1-c)*nx*ny	c+(1-c)*ny^2		-nx*s+(1-c)*nz*ny;
		-ny*s+(1-c)*nx*nz	nx*s+(1-c)*nz*ny	c+(1-c)*nz^2];
