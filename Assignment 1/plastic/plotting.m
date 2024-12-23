%figure(1); clear all; clc;load('plastic24micro.mat')
         subplot(1,3,1),cla
         plot_seff(node,elemnode,intpoints,seff);
         
         subplot(1,3,2),cla, mag=1e0;
         patch('vertices', node + mag*[un(1:2:end) un(2:2:end)],'faces', elemnode, 'edgecol', 'k', 'facecol', [.8,.9,1]); axis equal;
         axis([-1 1 -1 1]*L/2*1.1); axis equal; axis off; title(['deformation (magnified x', num2str(mag), ')'],'fontsize',14); drawnow
         
         subplot(1,3,3),cla,
         plot(tn*1e6,eff_stress*1e-6,'b-','linewidth',2); 
         k(1)=xlabel('t [µs]');k(2)=ylabel('\sigma_e [MPa]'); k(3)=title('Effective stress at crack edge'); set(k,'fontsize',14); grid on; drawnow      
            
%          i=B1;
%          subplot(2,2,1),cla,plot(tn*1e6,eff_stress_B1*1e-6,'b-','linewidth',2); 
%          k(1)=xlabel('t [µs]');k(2)=ylabel('\sigma_e [MPa]'); title(['\sigma_{eff} at B1 and t=' num2str(1e6*t) ' µs [MPa]']); set(k,'fontsize',14); grid on; drawnow
%          
%          i=B2;
%          subplot(2,2,2),cla,plot(tn*1e6,eff_stress_B2*1e-6,'b-','linewidth',2); 
%          k(1)=xlabel('t [µs]');k(2)=ylabel('\sigma_e [MPa]'); title(['\sigma_{eff} at B2 and t=' num2str(1e6*t) ' µs [MPa]']); set(k,'fontsize',14); grid on; drawnow
%          
%          i=B3;
%          subplot(2,2,3),cla,plot(tn*1e6,eff_stress_B3*1e-6,'b-','linewidth',2); 
%          k(1)=xlabel('t [µs]');k(2)=ylabel('\sigma_e [MPa]'); title(['\sigma_{eff} at B3 and t=' num2str(1e6*t) ' µs [MPa]']); set(k,'fontsize',14); grid on; drawnow
%          
%          i=B4;
%          subplot(2,2,4),cla,plot(tn*1e6,eff_stress_B4*1e-6,'b-','linewidth',2); 
%          k(1)=xlabel('t [µs]');k(2)=ylabel('\sigma_e [MPa]'); title(['\sigma_{eff} at B4 and t=' num2str(1e6*t) ' µs [MPa]']); set(k,'fontsize',14); grid on; drawnow
