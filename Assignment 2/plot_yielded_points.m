function plot_yielded_points(node,elemnode,intpoints,seffmax,s0)
% plot function

patch('vertices',node,'faces',elemnode,'facecolor','w','linestyle','-'); 
hold on, i=find(seffmax(:)>s0); x=intpoints(:,1,:); y=intpoints(:,2,:); plot(x(i),y(i),'r+'); 
title('yielded integration points','fontsize',14); axis equal; 
axis([min(node(:,1)) max(node(:,1)) min(node(:,2)) max(node(:,2))]); axis off; drawnow