function q_pos = MGI_pos(tbase_eff, est)

global axe1_d axe2_d axe3_d axe4_d


if nargin < 2
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

%tbase_eff
q_pos(1) = tbase_eff(3) - axe1_d + axe3_d;
q_pos(2) = tbase_eff(1); % - axe2_d;
q_pos(3) = tbase_eff(2) - axe2_d; % - axe3_d;