function [u,R] = solveSys(NdofsXnode,Ndofs,fixNod,KG,f)

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

u(vl) = ul;
u(vr) = ur;

R(vl) = 0;
R(vr) = Rr;
end