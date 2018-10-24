function f = computeF(NdofsXnode, Ndofs, Fext)
        Nforces = size(Fext, 1);
        f = zeros(Ndofs, 1);
        for k = 1:Nforces
            A = Fext(k,1);
            a = Fext(k,2);
            f(NdofsXnode*(A-1)+a,1) = f(NdofsXnode*(A-1)+a,1)+Fext(k,3);
        end
end