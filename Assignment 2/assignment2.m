% Simple FE code including J2-plasticity, bi-linear hardening
clc; clear, close all
         
%load a mesh
%[node,elemnode]=meshrect; 
[node,elemnode]=meshhole; 


% No of elements
Nel=size(elemnode,1); 
L=max(node(:,1)); %Length
W=max(node(:,2)); %Width

% Parameters
t=0.63e-3;    %thickness
du=0.2e-4; %incremental displacement load
u=0;       %set initial total displacement
f=0;       %set initial total force

% Material parameters
E=44e9; v=0.3; s0=190e6; n_exp=13; % E, v, yield stress, plastic modulus
% Sn youngs modulus
% Elastic constitutive matrix, plane stress
C = E/(1-v^2)*[1  v   0  
             v  1   0
             0  0 (1-v)/2];

% Define number of freedoms. Each node has 2 dofs
ndof=2*size(node,1);

% Allocate some variables/matrices
U=zeros(ndof,1);          % displacement of each dof (one node has 2 dofs: displacement in x-dir, displacement in y-dir)
stress=zeros(4,3,Nel);    % stress matrix [4 integration points, 3 components (sig_xx, sig_yy, tau_xy), in all Nel elements]
strain=zeros(4,3,Nel);    % strain matrix [4 integration points, 3 components (eps_xx, eps_yy, gam_xy), in all Nel elements]
seffmax=s0*ones(4,1,Nel); % present yield surface, initial and onset at s0
         
% Define global matrix of all integration points stiffness
% (note that each integration point has its own D!!)
Dep=repmat(C,[4 1 Nel]); %4 inegration points/element -> 4x(3x3) for each element

% Boundary nodes B1, B3 and their dofs b1x, b1y, b3x, b3y
%       +-------------B2-----------+
%       |                          |
%       B1                         B3
%       |                          |
%       +-------------B4-----------+
%                   
%   ^ y
%   |
%   +--> x
B1=find(node(:,1)==0); b1x=2*B1-1; b1y=2*B1;
B3=find(node(:,1)==L); b3x=2*B3-1; b3y=2*B3;

step=0; 
while u/L <= 0.032 %We load in xx steps

step=step+1;

% Boundary conditions
BC=[b1x    0*b1x
    b1y(1) 0
    b3x 0*b3x+du];

% Assemble the global stiffness matrix K
K=zeros(ndof,ndof);
for el=1:Nel
    %Compute material stiffness matrix Dep
    [Dep(:,:,el),seffmax(:,:,el)] = Constitutive_matrix_power_hardening(stress(:,:,el),seffmax(:,:,el),s0,E,v,n_exp, du); 
                                                                       
    n=elemnode(el,:);
    ex=node(n,1)'; ey=node(n,2)';
    Ke=plane_iso4kep(ex,ey,t,Dep(:,:,el));     
    elemdof=[n(1)*2-1 n(1)*2 n(2)*2-1 n(2)*2 n(3)*2-1 n(3)*2 n(4)*2-1 n(4)*2];   
    K(elemdof,elemdof)=K(elemdof,elemdof)+Ke;
end

% internal force vector
dF=zeros(ndof,1); %No internal tractions are given

% solve for incremental displacements!
[dU,dQ]=solveq(K,dF,BC);

% when increment is solved: extract incremental stresses, strains, integration points
for el=1:Nel
    n=elemnode(el,:);
    ex=node(n,1)'; ey=node(n,2)';
    elemdof=[n(1)*2-1 n(1)*2 n(2)*2-1 n(2)*2 n(3)*2-1 n(3)*2 n(4)*2-1 n(4)*2]; 
    ed=dU(elemdof)';
    [dstress(:,:,el),dstrain(:,:,el),intpoints(:,:,el)]=plane_iso4sp(ex,ey,Dep(:,:,el),ed);
end

% Update stresses and strains after load increment increment
stress=stress+dstress;
strain=strain+dstrain;

% Compute updated von Mises effective stress, if necessary. 2D plane stress!
seff=sqrt( stress(:,1,:).^2+stress(:,2,:).^2-stress(:,1,:).*stress(:,2,:)+3*stress(:,3,:).^2 );

% Update globale load and reaction forces
u(step+1)=u(step)+mean(dU(b3x)); % update displacement load
f(step+1)=f(step)+sum(dQ(b3x));  % update reaction force
U=U+dU;                          % update total displacements of all dofs

% draw some plots
figure(1)
subplot(2,2,1), cla, plot((u/L)*100,f/(t*W)*1e-6,'bo-'), grid on, xlabel('u/L','fontsize',14); ylabel('\sigma [MPa]','fontsize',14);
subplot(2,2,3), cla, plot_yielded_points(node,elemnode,intpoints,seffmax,s0) % /s0  /\sigma_0
magnificationfactor=4; 
subplot(2,2,4), cla, plot_deformed_mesh(node,elemnode,U,magnificationfactor)
subplot(2,2,2), cla, plot_plastic_zone(node,elemnode,intpoints,seff,s0)

end % return if u<L/50


