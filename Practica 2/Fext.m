function [Fext,acc] = fext(Nnodes, M, u, temps, Rhoa, Cd, g, m_p, m_r, Nelements, Tn)

%S'a de tenir en compte la força de les barres a cada node i la força del
%drag, la inercial i la gravitatoria. 
% Fext = [ 1              2              3
%         Node        Direccio        Magnitud    ]

 %Aqui definim l'inerval de temps i ens quedarem amb l'ultim t obtingut. 
Fext_g = zeros(Nnodes, 3, length(temps));
Fext_i = zeros(Nnodes, 3, length(temps));
Fext_d = zeros(Nnodes, 3, length(temps));
acc = zeros(length(temps));
for s = 1:length(temps);
    t=s;
    uc=vpa(subs(u));    
    D = 0.5*Rhoa*uc^2*Cd;
    a = g - D/(sum(M)+m_p+m_r);
    
    %Gravity force
for i=1:Nelements
    Fext_g(Tn(i,1),2,s) = 3; %Actuen tots en la direccio z
    Fext_g(Tn(i,2),2,s) = 3;
    Fext_g(Tn(i,1),3,s) = -M(i)*g/2;
    Fext_g(Tn(i,2),3,s) = -M(i)*g/2;
    Fext_g(Tn(i,1),1,s) = Tn(i,1);
    Fext_g(Tn(i,2),1,s) = Tn(i,2);
end

% Inertial Forces
for i=1:Nelements
    Fext_i(Tn(i,1),3,s) = -M(i)*a/2;
    Fext_i(Tn(i,2),3,s) = -M(i)*a/2;
end

% Drag Forces
F1 = D-m_r*g-m_r*a;
Fext_d(1,3,s) = -m_p*g;
Fext_d(6,3,s) = F1/16;
Fext_d(7,3,s) = F1/8;
Fext_d(8,3,s) = F1/16;
Fext_d(9,3,s) = F1/8;
Fext_d(10,3,s) = F1/4;
Fext_d(11,3,s) = F1/8;
Fext_d(12,3,s) = F1/16;
Fext_d(13,3,s) = F1/8;
Fext_d(14,3,s) = F1/16;
acc(s) = a;
end
Fext = Fext_g; %+ Fext_i + Fext_d;
end
