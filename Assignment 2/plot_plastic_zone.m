function plot_plastic_zone(node,elemnode,intpoints,seff,s0)

SE=seff(:); X=intpoints(:,1,:); X=X(:); Y=intpoints(:,2,:); Y=Y(:);

%interpolate value in integration points to nodal points
F = scatteredInterpolant(X,Y,SE,'natural'); se=F(node);

colormap jet
patch('vertices',node,'faces',elemnode,'facevertexcdata',se/s0,'facecolor','interp','linestyle','none');
caxis([-1 3]); title('\sigma_e/\sigma_0','fontsize',14); axis equal; 
axis([min(node(:,1)) max(node(:,1)) min(node(:,2)) max(node(:,2))]); axis off; colorbar; drawnow

