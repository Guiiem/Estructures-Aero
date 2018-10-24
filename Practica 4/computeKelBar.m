function Kel = computeKelBar(Nelements, NdofsXelement, x, Tn, Tmat, mat, E)

Kel = zeros(NdofsXelement, NdofsXelement, Nelements);

for e = 1:Nelements
    x1 = x(Tn(e,1),1);
    x2 = x(Tn(e,2),1);
    y1 = x(Tn(e,1),2);
    y2 = x(Tn(e,2),2);
    
    
    A = mat(Tmat(e),2);
    l = sqrt((x2-x1)^2+(y2-y1)^2);
    
    s = (y2-y1)/l; %sinus de l'angle que forma la barra
    c = (x2-x1)/l; %cosinus de l'angle que forma la barra
    
    
    a = (6*mat(Tmat(e),1))/l;
    b = (12*mat(Tmat(e),1))/(l^2);
    
    
    Kel(:,:,e) = (E/l)*[A*(c^2)+b*(s^2) (A-b)*c*s -a*s -(A*(c^2)+b*(s^2)) -(A-b)*c*s -a*s;
        (A-b)*c*s A*(s^2)+b*(c^2) a*c -(A-b)*c*s -(A*(s^2)+b*(c^2)) a*c;
        -a*s a*c 4*mat(Tmat(e),1) a*s -a*c 2*mat(Tmat(e),1);
        -(A*(c^2)+b*(s^2)) -(A-b)*c*s a*s A*(c^2)+b*(s^2) (A-b)*c*s a*s;
        -(A-b)*c*s -(A*(s^2)+b*(c^2)) -a*c (A-b)*c*s A*(s^2)+b*(c^2) -a*c;
        -a*s a*c 2*mat(Tmat(e),1) a*s -a*c 4*mat(Tmat(e),1)];
end
end