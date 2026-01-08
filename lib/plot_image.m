function plot_image(curr_image, camera, col)


plot(curr_image(1,:), curr_image(2,:), col);

axis([camera.corner_im(1) camera.corner_im(2) camera.corner_im(3) camera.corner_im(4)]);

drawnow;