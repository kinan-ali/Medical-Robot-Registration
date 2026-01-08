function [im_pts, liste_visible_points] = camera_proj(Tcam_obj, camera, object, noise_ampl)

fprintf(1, 'Image computation\n');

if nargin < 4
	noise_ampl = 0;
end

%object

nb_pts_obj = size(object, 2);

liste_visible_points = zeros(size(object,2), 1);

Rcam_obj = Tcam_obj(1:3, 1:3);
tcam_obj = Tcam_obj(1:3, 4);

P_in_cam = Rcam_obj * object + repmat(tcam_obj, 1, size(object,2));


% Points projetes sur le capteur
m = [P_in_cam(1,:)./P_in_cam(3,:); P_in_cam(2,:)./P_in_cam(3,:); ones(1,size(object,2))];

% points image obtenus
im_pts = camera.K*m + [2*noise_ampl*(rand(2,nb_pts_obj)-0.5*ones(2,nb_pts_obj)); zeros(1,nb_pts_obj)];


% Test si les points sont dans l'image
for kk=1:nb_pts_obj
	if((test_pix([im_pts(1,kk); im_pts(2,kk)], camera.corner_im(1), camera.corner_im(2), camera.corner_im(3), camera.corner_im(4)) > 0) && (P_in_cam(3,kk) > 0))
		liste_visible_points(kk) = 1;
	else
		liste_visible_points(kk) = 0;
		im_pts(:,kk) = [-1000; -1000; -1000];
	end


end

fprintf(1, 'Image computed\n');
