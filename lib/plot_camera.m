function plot_camera(not_used, Tworld_cam, col, flag_fustrum, flag_count)

persistent center_prev;

if nargin < 4
	flag_count = 1;
	if nargin < 3
		flag_fustrum = 0;
		if nargin < 2
			col = 'k'	
		end	
	end
end

ang_view = pi/2;
scale = 50;
corner = scale*sin(ang_view/2);

camera_origin = [[corner; -corner; scale], [corner; corner; scale], [-corner; corner;scale], [-corner; -corner;scale], [corner; -corner; scale]]; 

center = Tworld_cam(1:3,4);
Rworld_cam = Tworld_cam(1:3, 1:3);

camera = Rworld_cam*camera_origin + repmat(center, 1, 5);

if flag_fustrum==1
	for kk=1:4
		plot3([camera(1,kk), camera(1,kk+1)],[camera(2,kk), camera(2,kk+1)],[camera(3,kk), camera(3,kk+1)], col); 
		hold on;
	end
	
	for kk=1:4
		plot3([center(1), camera(1,kk)],[center(2), camera(2,kk)],[center(3), camera(3,kk)], col); 
		hold on;
	end
end

xcross = sprintf('%s+', col);
lline = sprintf('%s-', col);
if ((flag_count == 1) || isempty(center_prev))
	P = plot3(center(1), center(2),center(3), xcross);
else
	P = plot3([center_prev(1), center(1)], [center_prev(2),center(2)],[center_prev(3), center(3)], lline);
end

if flag_count == -1 % fin de trajectoire
	center_prev = [];
else
	center_prev = center;
end

set(P, 'linewidth', 2);

