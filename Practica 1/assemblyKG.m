function KG = assemblyKG(Nelements,NdofsXelement,Ndofs,Td,Kel)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  Nelements       Total number of elements
%                  NdofsXelement   Number of DOFs per element
%                  Ndofs           Total number of DOFs
%   - Td    DOFs connectivities table [Nelements x NdofsXelement]
%            Td(e,i) - DOF i associated to element e
%   - Kel   Elemental stiffness matrices [NdofsXelement x NdofsXelement x Nelements]
%            Kel(i,j,e) - Term in (i,j) position of stiffness matrix for element e
%--------------------------------------------------------------------------
% It must provide as output:
%   - KG    Global stiffness matrix [Ndofs x Ndofs]
%            KG(I,J) - Term in (I,J) position of global stiffness matrix
%--------------------------------------------------------------------------
KG = zeros(Ndofs,Ndofs);

for e=1:Nelements
    for i=1:NdofsXelement
        I=Td(e,i);
        for j=1:NdofsXelement
            J=Td(e,j);
            KG(I,J)=KG(I,J)+Kel(i,j,e);
        end
    end
end