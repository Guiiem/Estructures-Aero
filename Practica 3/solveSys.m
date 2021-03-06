function [u,R] = solveSys(NdofsXnode,Ndofs,fixNod,KG,f)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  NdofsXnode      Number of DOFs per node
%                  Ndofs           Total number of DOFs
%   - fixNod  Prescribed displacements data [Npresc x 3]
%              fixNod(k,1) - Node at which the some DOF is prescribed
%              fixNod(k,2) - DOF (direction) at which the prescription is applied
%              fixNod(k,3) - Prescribed displacement magnitude in the corresponding DOF
%   - KG      Global stiffness matrix [Ndofs x Ndofs]
%              KG(I,J) - Term in (I,J) position of global stiffness matrix
%   - f       Global force vector [Ndofs x 1]
%              f(I) - Total external force acting on DOF I
%--------------------------------------------------------------------------
% It must provide as output:
%   - u       Global displacement vector [Ndofs x 1]
%              u(I) - Total displacement on global DOF I
%   - R       Global reactions vector [Ndofs x 1]
%              R(I) - Total reaction acting on global DOF I
%--------------------------------------------------------------------------
% Hint: Use the relation between the DOFs numbering and nodal numbering to
% determine at which DOF in the global system each displacement is prescribed.
Npresc=size(fixNod);
Npresc=Npresc(1,1);

vr = zeros(1,Npresc);
ur = zeros(Npresc,1);

for k=1:Npresc
    A=fixNod(k,1);
    i=fixNod(k,2);
    vr(1,k)=NdofsXnode*(A-1)+i;
    ur(k,1)=fixNod(k,3);
end

vl=setdiff(1:Ndofs,vr);

Kll = KG(vl,vl);
Klr = KG(vl,vr);
Krl = KG(vr,vl);
Krr = KG(vr,vr);

fl = f(vl,1);
fr = f(vr,1);

ul = Kll\(fl - Klr*ur);

Rr = Krr*ur + Krl*ul - fr;

u(vl,1) = ul;
u(vr,1) = ur;

R(vl,1) = 0;
R(vr,1) = Rr;
end
