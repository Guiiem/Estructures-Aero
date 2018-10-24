function Kel = computeKelBar(Nelements, NdofsXelement, x, Tn, Tmat, mat, E)
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
Kel = zeros(NdofsXelement, NdofsXelement, Nelements);
for e = 1:Nelements
    x1 = x(Tn(e,1),1);
    x2 = x(Tn(e,2),1);
    y1 = x(Tn(e,1),2);
    y2 = x(Tn(e,2),2);
    
    l = sqrt((x2-x1)^2+(y2-y1)^2);
    
    s = (y2-y1)/l;
    c = (x2-x1)/l;
    
    alpha = (6*mat(Tmat(e),1))/l;
    beta = (12*mat(Tmat(e),1))/(l^2);
    A = mat(Tmat(e),2);
    
    Kel(:,:,e) = (E/l)*[A*(c^2)+beta*(s^2) (A-beta)*c*s -alpha*s -(A*(c^2)+beta*(s^2)) -(A-beta)*c*s -alpha*s
        (A-beta)*c*s A*(s^2)+beta*(c^2) alpha*c -(A-beta)*c*s -(A*(s^2)+beta*(c^2)) alpha*c
        -alpha*s alpha*c 4*mat(Tmat(e),1) alpha*s -alpha*c 2*mat(Tmat(e),1)
        -(A*(c^2)+beta*(s^2)) -(A-beta)*c*s alpha*s A*(c^2)+beta*(s^2) (A-beta)*c*s alpha*s
        -(A-beta)*c*s -(A*(s^2)+beta*(c^2)) -alpha*c (A-beta)*c*s A*(s^2)+beta*(c^2) -alpha*c
        -alpha*s alpha*c 2*mat(Tmat(e),1) alpha*s -alpha*c 4*mat(Tmat(e),1)];
end
end