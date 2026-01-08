
function plot_patient_body()

global points px py pz pc

figure_handle = gcf;

% comment/uncomment the following lines as you wish:

% plot 3D points
  %plot3(points(:,1),points(:,2),points(:,3),'bx','markersize',1)
    
  % plot polygon's lines 
  %  line(px(:,:)',py(:,:)',pz(:,:)','color','b');
    
   % fills a polygon with a given color
   pa = patch(1.5*px(:,:)',(-1.5*py(:,:)')-200,2.5*pz(:,:)',pc(:,:)');
   %pa = patch(px(:,:)',py(:,:)',pz(:,:)',pc(:,:)');
   % If using this, you might want to render in z-buffer
   %  mode for lighting effects
   set(figure_handle,'Renderer','zbuffer')
   light('color','w')
   %light('color','c','position',[-1 5 0.2])

    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    
%shading flat                                                         % select one for use with 
shading interp                                                       % the patch option

axis equal
    
pa.FaceAlpha = 0.4;    % Set constant transparency 

%Original nancy in VRML by Cindy Ballreich (cindy@ballreich.net) 3NAME3D.
%Matlab conversion by Ben Tordoff and Walterio Mayol RRG Univ. of Oxford.
%see copyright notice at the top of file.
