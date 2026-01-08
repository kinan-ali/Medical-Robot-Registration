function DisplayImage

global im_pts im_pts_mark_inst camera_real

title('position of target in current image');
plot_image(im_pts, camera_real, 'r*')
hold on;

if(~isempty(im_pts_mark_inst))
    plot_image(im_pts_mark_inst(:,1), camera_real, 'g*')
    legend('target', 'instrument tip');
else
  legend('target');
end
