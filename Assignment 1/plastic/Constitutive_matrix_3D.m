function [D,seffmax]=Constitutive_matrix_3D(stress,seffmax,s0,E,v,H)

%Define 6x6 (3D) elastic constitutive matrix
C=E/(1+v)/(1-2*v)*[1-v  v    v    0           0       0
                    v  1-v   v    0           0       0
                    v   v   1-v   0           0       0
                    0   0    0 (1-2*v)/2      0       0
                    0   0    0    0      (1-2*v)/2    0
                    0   0    0    0           0   (1-2*v)/2];
                
D=[];
for in=1:4 %loop over all 4 integration points     
    sig11=stress(in,1);
    sig22=stress(in,2);
    sig12=stress(in,3);
    sig33=0; % PLANE STRESS
    sig13=0; % PLANE STRESS
    sig23=0; % PLANE STRESS    
    sh=(sig11+sig22+sig33)/3; %hydrostatic stress    
    
    %deviator stresses sij
    s11=sig11-sh;
    s22=sig22-sh;
    s33=sig33-sh;
    s12=sig12;
    s21=sig12;
    s23=sig23;
    s13=sig13;
    
    %J2 effective stress, plane case
    se=sqrt(3/2*(s11.*s11 + s22.*s22 + s33.*s33 + 2*s12.*s12));      

    if se >= seffmax(in) % yielding 
           beta = 1;
    else                 % no yielding
           beta = 0; 
    end

    seffmax(in) = max([seffmax(in) se]); %Update seffmax, contains highest value of se ever reached
    Dp=zeros(6,6);      % 3-D plastic stiffness matrix
    
    % Compute the plastic part of the constitutive matrix
    if beta==1  
        Et=H*E/(H+E); E_Et=E/Et;  % bi-linear plasticity, see lecture         
        h=9/(4*E*se^2)*(E_Et-1);
        sij=[s11 s22 s33 s12 s23 s13]'; %stress deviator vector

        Dp=h*C*sij*sij'*C/(1+h*sij'*C*sij);        
    
    end % end plastic part (i.e if beta=1)
    Dep=C-Dp; % elastic-plastic incremental tangential stiffness matrix 3D

    % 3-D --> 2-D (we have PLANE STRESS, but all computations are 3D)
    % format is [11 22 33 12 23 13] only 11, 22, 12 are needed

    invDep=inv(Dep); Dep2D=inv( invDep([1 2 4],[1 2 4]) ); 

    D=[D; Dep2D]; %assemble big D=[D1; D2; ....; D(Nin*Nin)]   
end
     

    
