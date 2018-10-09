function [eps,sig,eps_t,sig_t] = computeStrainStressBar(Ndim,Nelements,u,Td,x,Tn,mat,Tmat,temps)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  Ndim        Problem's dimensions
%                  Nelements   Total number of elements
%   - u     Global displacement vector [Ndofs x 1]
%            u(I) - Total displacement on global DOF I
%   - Td    DOFs connectivities table [Nelements x NdofsXelement]
%            Td(e,i) - DOF i associated to element e
%   - x     Nodal coordinates matrix [Nnodes x Ndim]
%            x(a,i) - Coordinates of node a in the i dimension
%   - Tn    Nodal connectivities table [Nelements x NnodesXelement]
%            Tn(e,a) - Nodal number associated to node a of element e
%   - mat   Material properties table [Nmat x NpropertiesXmat]
%            mat(m,1) - Young modulus of material m
%            mat(m,2) - Section area of material m
%   - Tmat  Material connectivities table [Nelements]
%            Tmat(e) - Material index of element e
%--------------------------------------------------------------------------
% It must provide as output:
%   - eps   Strain vector [Nelements x 1]
%            eps(e) - Strain of bar e
%   - sig   Stress vector [Nelements x 1]
%            sig(e) - Stress of bar e
%--------------------------------------------------------------------------
eps_t  = zeros(Nelements,1,length(temps));
sig_t = zeros(Nelements,1,length(temps));
eps = zeros(Nelements,1);
sig = zeros(Nelements,1);
de=zeros(2*Ndim,1);
for s = 1:length(temps)
for e=1:Nelements
    x1e=x(Tn(e,1),1);
    y1e=x(Tn(e,1),2);
    z1e=x(Tn(e,1),3);
    x2e=x(Tn(e,2),1);
    y2e=x(Tn(e,2),2);
    z2e=x(Tn(e,2),3);
    
    le=sqrt((x1e-x2e)^2+(y1e-y2e)^2+(z2e-z1e)^2);
    Re=1/le*[
       x2e-x1e, y2e-y1e,z2e-z1e,0,0,0;
       0,0,0,x2e-x1e, y2e-y1e,z2e-z1e;
       ];
       
    %Tot el qu ehem calculat anteriorment no depen del temps: le i Re. Re
    %[2x6]
%obtain displacements
for r=1:6
       p=Td(e,r);
       de(r)=u(p,1,s); 
end
dprim=Re*de; %ho calculem per al temps s.
%calculate Strain and Stress
eps(e,1)=(1/le)*[-1 1]*dprim; 
sig(e,1)=mat(Tmat(e),1)*eps(e,1);
end
for index_1 = 1:Nelements
eps_t(index_1,1,s) = eps(index_1,1); %Ens guardem els valors per posterior
sig_t(index_1,1,s) = sig(index_1,1); % ús en una matriu.
end
end
end
