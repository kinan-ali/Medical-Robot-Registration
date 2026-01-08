function Trob_eff = MGD(q, tbase_troc, est)

global axe1_d axe2_d axe3_d axe4_d

%tbase_troc

if nargin < 3
	% Param fixes exacts
	a = 500; %100;
	b = 0; % 200
	c = 0; %200;
	d = -pi/2; %pi/24;
	e = 0; %pi/24;
	f = 0; %pi/24;

else
	% Param fixes estimes
	a = 100 *0.9; %100;
	b = 10; % 200
	c = 10; %200;
	d = -pi/2 + pi/8; %pi/24;
	e = pi/8; %pi/24;
	f = pi/8; %pi/24;

	
	%a = 110;
	%b = 210;
	%c = 190;
	%d = 0;
	%e = 0;
	%f = 0;

	%a = 100;
	%b = 200;
	%c = 202;
	%d = pi/24;
	%e = pi/24;
	%f = 0;


end


% MGD d'un robot cartesien + poignet

%axe1_d
%q

H01 = [eye(3) [0; 0; q(1) + axe1_d];
	0 0 0 1];

H12 = [theta2r([0; pi/2;0]) [q(2); axe2_d; 0];
	0 0 0 1];

%H02 = H01*H12

	
H23 = [theta2r([-pi/2;0;0]) [axe3_d; q(3); 0]
	0 0 0 1];	 
	
%H03 = H02*H23	
	
H34 = [theta2r([0;pi/2; 0;0]) [0; axe4_d; 0];
	0 0 0 1];	 
		
		
H04 = H01 * H12*H23*H34;		
		
tbase_eff = H04(1:3,4);


eff_troc_in_base = tbase_troc - tbase_eff;
axe = eff_troc_in_base / norm(eff_troc_in_base);

% 0;0;-1 est l'orientation naturelle de l'instrument sans contrainte de trocart
sang_u = cross([0;0;-1], axe);
sang = norm(sang_u);

if axe(3) <= 0
	ang = asin(sang);
else
	ang = pi - asin(sang);
end

if abs(sang) > 1e-5
	u = sang_u / sang;
else
	u = [0;0;1]; 
end

Rbase_eff_nom = theta2r([pi;0;0]);		

u_in_eff_nom = Rbase_eff_nom' * u;
Reff_nom_eff = thetau2r(ang*u_in_eff_nom);		
				
Rbase_eff = Rbase_eff_nom*Reff_nom_eff;		

Trob_eff = [Rbase_eff, tbase_eff; 0 0 0 1];
		
%H47 = [theta2r([q(4)+d 0 0])*theta2r([0 q(5)+e 0])* theta2r([0 0 q(6)+f]) zeros(3,1);
%	0 0 0 1]; 

%Trob_eff = H12*H23*H34*H47;
