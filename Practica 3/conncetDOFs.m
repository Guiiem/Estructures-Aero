function Td = connectDOFs(Nelements,NnodesXelement,NdofsXnode,Tn)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  Nelements       Total number of elements
%                  NnodesXelement  Number of nodes per element
%                  NdofsXnode      Number of DOFs per node
%   - Tn    Nodal connectivities table [Nelements x NnodesXelement]
%            Tn(e,a) - Nodal number associated to node a of element e
%--------------------------------------------------------------------------
% It must provide as output:
%   - Td    DOFs connectivities table [Nelements x NdofsXelement]
%            Td(e,i) - DOF i associated to element e
%--------------------------------------------------------------------------
% Hint: Use the relation between the DOFs numbering and nodal numbering.
Td = zeros(Nelements, NnodesXelement*NdofsXnode);
for i = 1:Nelements
    for j = 1:NnodesXelement
       for k = 1:NdofsXnode
        Td(i,j*NdofsXnode+k-NdofsXnode) = Tn(i,j)*NdofsXnode+k-NdofsXnode;
       end
    end
end


end