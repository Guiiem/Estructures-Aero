function Td = connectDOFs(Nelements, NnodesXelement, NdofsXnode, Tn)
Td = zeros(Nelements, NnodesXelement*NdofsXnode);

for e = 1:Nelements
    for j = 1:NnodesXelement
        for k = 1:NdofsXnode
        	Td(e, (j-1)*NdofsXnode+k) = Tn(e,j)*NdofsXnode+k-NdofsXnode;
        end
    end
end
end