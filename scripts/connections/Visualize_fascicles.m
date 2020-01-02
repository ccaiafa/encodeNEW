% Below are a set of local matlab functions that are sued in this script.
function [] = Visualize_fascicles(fe, fascicles_ind, voxelA_ind, voxelB_ind, fig_name)
% 
% This function is used to visualize the anatomy of part of connectome
% fascicles and a region of interest (ROI). 
% 
% It calls functions from github.com/francopestilli/mba

colors     = {[.1 .25 .65]};
viewCoords = [0,0];

fg{1}          = feGet(fe,'fibers img');
fg{1}.fibers   = fg{1}.fibers(fascicles_ind);

% plot fascicles
[~, ~] = plotFasciclesNoAnat(fg, colors, viewCoords, fig_name, [1]);

set(gcf,'Color',[1 1 1])
hold on
offset = 1.5;
% Plot region of interest (ROI), anatomy voxels region A.
plot3(fe.roi.coords(voxelA_ind,1)-offset, ...
         fe.roi.coords(voxelA_ind,2)-offset, ...
         fe.roi.coords(voxelA_ind,3)-offset, 'ro', ...
         'markerfacecolor','r', 'markersize',10)
% Plot region of interest (ROI), anatomy voxels region B.
plot3(fe.roi.coords(voxelB_ind,1)-offset, ...
         fe.roi.coords(voxelB_ind,2)-offset, ...
         fe.roi.coords(voxelB_ind,3)-offset, 'bo', ...
         'markerfacecolor','b', 'markersize',10)

end

function [fig_h, light_h] = plotFasciclesNoAnat(fascicles, color, viewCoords, fig_name,tracts_to_clean)
% 
% This function is used to visualize the anatomy of part of connectome
% fascicles. 
%
% It calls functions from github.com/francopestilli/mba

fig_h = figure('name',fig_name,'color','k');
hold on
set(gca,'visible','off','Color','w')
for iFas  = 1:length(tracts_to_clean)
    [~, light_h] = mbaDisplayConnectome(fascicles{ tracts_to_clean(iFas) }.fibers,fig_h,color{ tracts_to_clean(iFas) },'single');
    delete(light_h)
end
view(viewCoords(1),viewCoords(2))
light_h = camlight('right');
lighting phong;
set(gcf,'Color',[1 1 1])
drawnow

end
