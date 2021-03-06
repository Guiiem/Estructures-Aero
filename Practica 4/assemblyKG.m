function KG = assemblyKG(Ndofs, NdofsXelement, Nelements, Td, Kel)

KG = zeros(Ndofs,Ndofs);

for e=1:Nelements
    for i=1:NdofsXelement
        I=Td(e,i);
        for j=1:NdofsXelement
            J=Td(e,j);
            KG(I,J)=KG(I,J)+Kel(i,j);
        end
    end
end

end