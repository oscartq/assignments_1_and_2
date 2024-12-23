function plot_deformed_mesh(node,elemnode,U,factor)

Ux=U(1:2:end); Uy=U(2:2:end); %Nodal displacement

X=[node+factor*[Ux Uy]];
patch('vertices',X,'faces',elemnode,'facecolor','b','linestyle','-','facealpha',0.3); 
title(['deformed (x' num2str(factor) ') mesh'],'fontsize',14); axis equal;
axis([min(X(:,1)) max(X(:,1)) min(X(:,2)) max(X(:,2))]); axis off;  drawnow