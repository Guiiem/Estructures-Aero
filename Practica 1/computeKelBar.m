function [Kel,massa_bc] = computeKelBar(Ndim,Nelements,x,Tn,mat,Tmat)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  Ndim        Problem's dimensions
%                  Nelements   Total number of elements
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
%   - Kel   Elemental stiffness matrices [NdofsXelement x NdofsXelement x Nelements]
%            Kel(i,j,e) - Term in (i,j) position of stiffness matrix for element e
%--------------------------------------------------------------------------
Kel = zeros(2*Ndim, 2*Ndim, Nelements);
massa_bc = 0;

for e = 1:Nelements
    x_1 = x(Tn(e,1),1);
    x_2 = x(Tn(e,2),1);
    y_1 = x(Tn(e,1),2); 
    y_2 = x(Tn(e,2),2);
    %%if Ndim == 3 %Per pode canviar-ho com volguem.
    z_1 = x(Tn(e,1),3);
    z_2 = x(Tn(e,2),3);
    %%else
    
    l = sqrt((x_2-x_1)^2+(y_2-y_1)^2+(z_2-z_1)^2);
    R_2= [x_2-x_1 y_2-y_1 z_2-z_1 0 0 0; 
        0 0 0 x_2-x_1 y_2-y_1 z_2-z_1];
    R = (1/l)*R_2; %2x6 matrix
    Kt = ((mat(Tmat(e),2)*mat(Tmat(e),1))/l)*[1, -1;-1, 1]; %2x2
    K = transpose(R)*Kt*R; %6x2 * 2x2 * 2x6; 
    massa_bc = massa_bc + mat(Tmat(e),2)*mat(Tmat(e),3)*l;
    for i = 1:4
        for j = 1:4
            Kel(i,j,e) = K(i,j);
        end
    end
end


end