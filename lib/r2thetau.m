function [theta,u] = r2thetau(R)

% Use : [theta,u] = r2thetau(R)
%
% Cette fonction exprime une rotation donnee sous
% forme matricielle sous la forme angle, axe
% Attention le resultat obtenu  n'est pas la seule solution possible
% (c.f. Vision par ordinateur p.240)
% On choisit ici d'avoir theta dans [0 +pi]

%fprintf(1, 'r2thetau');

temp = 0.5*(R(1,1) + R(2,2) + R(3,3) - 1);
if temp >= 1-1e-8
	temp = 1;
	theta = 0;
	u = [1;0;0];
	return;
end
if temp <= -1+1e-8
	temp = -1;
	theta = pi;
	[V,D] = eig(R)
	[r,c,v] = find(D > 0.5)
	u = V(1:3,r)
	return;	
end

theta = acos(temp);

%if(abs(theta) < 1e-8) 
%	theta = 0;
%	u = [1; 0; 0];
%else

sintheta = sin(theta);

n1 = (R(3,2) - R(2,3)) / (2*sintheta);
n2 = -(R(3,1) - R(1,3))/(2*sintheta);
n3 = (R(2,1) - R(1,2)) / (2*sintheta);

u = [n1;n2;n3];
%end
