function f = computeF(NdofsXnode, Ndofs, Fext)
        Nforces = size(Fext, 1);
        f = zeros(Ndofs, 1);
        for k = 1:Nforces
            node = Fext(k,1);
            dof = Fext(k,2);
            f(NdofsXnode*(node-1)+dof,1) = f(NdofsXnode*(node-1)+dof,1)+Fext(k,3);
        end
end