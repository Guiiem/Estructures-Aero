function f = computeF(NdofsXnode,Ndofs,Fext)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  NdofsXnode      Number of DOFs per node
%                  Ndofs           Total number of DOFs
%   - Fext  External nodal forces [Nforces x 3]
%            Fext(k,1) - Node at which the force is applied
%            Fext(k,2) - DOF (direction) at which the force acts
%            Fext(k,3) - Force magnitude in the corresponding DOF
%--------------------------------------------------------------------------
% It must provide as output:
%   - f     Global force vector [Ndofs x 1]
%            f(I) - Total external force acting on DOF I
%--------------------------------------------------------------------------
% Hint: Use the relation between the DOFs numbering and nodal numbering to
% determine at which DOF in the global system each force is applied.
f = zeros(Ndofs, 1);
for i=1:size(Fext,1);
if Fext(i,3)~=0
    if Fext(i,2) == 1
        f(Fext(i,1)*2-1)=Fext(i,3);
    else
        f(Fext(i,1)*2)=Fext(i,3);
end
end
end
