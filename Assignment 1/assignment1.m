% Assignment 1
clc; clear, close all

% Parameters
% Dimensions
L = 100e-3; % Length
H = 20e-3; % Height divded by two
thickness = 2e-3;   % Thickness

% Parameters
v   = 0.2; % Poisson's ratio
rho = 2450; % Density
E   = 32e9; % Young's modulus
c   = 15e6; % Viscous damping
c0  = sqrt(E/rho); % Longitudinal wave speed
t   = 0; % Time
dt  = 2e-8; % Time increment
F = 1e5; % Initial force

% Constitutive matrix
D = E/(1-v^2)*[1  v   0  %Plane stress
               v  1   0
               0  0 (1-v)/2];

% Get the square mesh
% nodes(x,y) and topology

% Geometry
%
%
%        B2
%     +-------+
% B1  |       | B3
%     +-------+
%        B4
%          

%load mesh
NelL=100;
[node,elemnode,B1,B2,B3,B4] = get_mesh_filter(L,H,NelL);

node(:,2) = node(:,2) + H/2;

B0 = find(node(:,1) >= -1e-10 & node(:,2) < 1e-10); %Take everything near 0 in height and at half length
%B2 = find(node(:,2) > H-1e-10); %Take all nodes that are close to the height
b0y = 2*B0; b0x = 2*B0-1; 
b2y = 2*B2; b2x = 2*B2-1;

BC = [b0y b0y*0];

% Define number of freedoms. Here each node have 2 dofs
ndof =  2*size(node,1);  % Total degrees of freedom
Nel  =  size(elemnode,1); % Number of elements
un   =  zeros(ndof,1);   % Displacements
vn   =  zeros(ndof,1);   % Velocities
an   =  zeros(ndof,1);   % Accelerations
strain = zeros(4,3,Nel);   % Strain
stress = zeros(4,3,Nel);   % Stress
intpoints = zeros(4,2,Nel);% Integration points


% Define global load vectors Fb and Fl
Fb = zeros(ndof,1);
Fl = zeros(ndof,1);

dL = L/(length(B2)-1);
Fb(b2y) = F*dL; Fb(b2y(1))= F*dL/2; Fb(b2y(end))= F*dL/2; 


% Mass matrix M and dampning matrix C
% Use advantage of exact same size of elements
% area Ae of the elements (equal!)
Ae=(L/(length(B2)-1))*(H/(length(B1)-1)); Mi=zeros(size(node,1),1); 

for i=1:size(node,1)
    j=find(elemnode(:)==i); 
    Mi(i)=rho*thickness*length(j)*Ae/4;
end
% Global mass matrix
M = sparse(1:ndof); M(1:2:end)=Mi; M(2:2:end)=Mi; M=diag(M);

% Global damping matrix
C = M * c/rho;

step = 0; 
while 1

    step = step + 1; 
    t = t + dt;
    
    %1. update displacements
    an(b0y) = 0;
    vn(b0y) = 0;
    
    un = un + dt*vn + 0.5*dt^2*an; % newmark eq 7
    

    %2. update internal forces and stiffness matrix 
    Fint = zeros(ndof,1);
    
    for el=1:Nel 
       
        m = elemnode(el,:); 
        n = [2*m(1)-1 2*m(1) 2*m(2)-1 2*m(2) 2*m(3)-1 2*m(3) 2*m(4)-1 2*m(4)];
        [Fe, Ke, strain(:,:,el), stress(:,:,el), intpoints(:,:,el)] = plane_4_iso(node(m,1)', node(m,2)', D, un(n)', thickness);
        Fint(n) = Fint(n) + Fe;
    end  
    
    %3. update velocities and accelerations: fully explicit method
    f = Fb - Fint - C*(vn + 0.5*dt*an); % newmark eq 8
    
    % New values for displacement
    an1 = (M + 0.5*dt*C) \ f; % newmark eq 9 but positive?
   
    vn = vn + dt*0.5*an + dt*0.5*an1; % newmark eq 10
    an = an1; 
    
    % Calculate effective stress
    seff = sqrt(stress(:,1,:).^2 + stress(:,2,:).^2 - stress(:,1,:).*stress(:,2,:) + 3*stress(:,3,:).^2 ); % calculateeffective stress 2D

    %4. Save some selected variables
    p(step+1) = Fint(b0y(1)) / (L/NelL*thickness);
    u(step+1) = un(b0y(1));
    tn(step+1) = t;
    dunorm(step+1) = norm(dt*vn + 0.5*dt^2*an)/norm(un);
    
    [row,col] = find(elemnode == B0(1)); % Where is the crack edge?
    eff_stress(step+1) = seff(col(1),:,row(1));

    %boundaries
    [sx,sy,sxy]=interpolate2nodes(stress,intpoints,node);
    seff_node = sqrt(sx.^2 + sy.^2 - sx.*sy + 3*sxy.^2);
    
    i = B1; eff_stress_B1(step+1) = mean(seff_node(i));
    i = B2; eff_stress_B2(step+1) = mean(seff_node(i));
    i = B3; eff_stress_B3(step+1) = mean(seff_node(i));
    i = B4; eff_stress_B4(step+1) = mean(seff_node(i));
    %%
    plot_step=100; % plot at every x step
      if rem(step,plot_step)==0 
         if step == plot_step, fprintf(1, '\nSteps done: '); end
         if rem(step,10*plot_step) == 0,fprintf(1,'\n '); end
         fprintf(1,'%1.0f, ',step);

         figure(1);
         subplot(2,4,1),cla
         plot_seff(node,elemnode,intpoints,seff);
         
         subplot(2,4,2),cla, mag=1e1;
         patch('vertices', node + mag*[un(1:2:end) un(2:2:end)],'faces', elemnode, 'edgecol', 'k', 'facecol', [.8,.9,1]); axis equal;
         axis([-1 1 -1 1]*L/2*1.1); axis equal; axis off; title(['deformation (magnified x', num2str(mag), ')'],'fontsize',14); drawnow
         
         subplot(2,4,3),cla,
         plot(tn*1e6,eff_stress*1e-6,'b-','linewidth',2); 
         k(1)=xlabel('t [µs]');k(2)=ylabel('\sigma_e [MPa]'); k(3)=title('Effective stress at crack edge'); set(k,'fontsize',14); grid on; drawnow      
            
         i=B1;subplot(2,4,5),cla,plot(tn*1e6,eff_stress_B1*1e-6,'b-','linewidth',2); 
         k(1)=xlabel('t [µs]');k(2)=ylabel('\sigma_e [MPa]'); title(['\sigma_{eff} at B1 and t=' num2str(1e6*t) ' µs [MPa]']); set(k,'fontsize',14); grid on; drawnow
         
         i=B2;subplot(2,4,6),cla,plot(tn*1e6,eff_stress_B2*1e-6,'b-','linewidth',2); 
         k(1)=xlabel('t [µs]');k(2)=ylabel('\sigma_e [MPa]'); title(['\sigma_{eff} at B2 and t=' num2str(1e6*t) ' µs [MPa]']); set(k,'fontsize',14); grid on; drawnow
         
         i=B3;subplot(2,4,7),cla,plot(tn*1e6,eff_stress_B3*1e-6,'b-','linewidth',2); 
         k(1)=xlabel('t [µs]');k(2)=ylabel('\sigma_e [MPa]'); title(['\sigma_{eff} at B3 and t=' num2str(1e6*t) ' µs [MPa]']); set(k,'fontsize',14); grid on; drawnow
         
         i=B4;subplot(2,4,8),cla,plot(tn*1e6,eff_stress_B4*1e-6,'b-','linewidth',2); 
         k(1)=xlabel('t [µs]');k(2)=ylabel('\sigma_e [MPa]'); title(['\sigma_{eff} at B4 and t=' num2str(1e6*t) ' µs [MPa]']); set(k,'fontsize',14); grid on; drawnow

      end

     %if dunorm(end) < 1e-5, break, end
     if t >= 1500e-6, break, end
end

fprintf(1,'\n\tFINISHED. TOTAL %1.0f steps made. \n',step);
