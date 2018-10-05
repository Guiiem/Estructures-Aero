function [Fext, acc] = fext(Nnodes, M, u, temps, Rhoa, Cd, g, m_p, m_r, Nelements, Tn, m_t, S)

%S'a de tenir en compte la for�a de les barres a cada node i la for�a del
%drag, la inercial i la gravitatoria. 
% Fext = [ 1              2              3
%         Node        Direccio        Magnitud    ]

 %Aqui definim l'inerval de temps i ens quedarem amb l'ultim t obtingut. 
Fext_g = zeros(Nnodes, 3, length(temps));
Fext_i = zeros(Nnodes, 3, length(temps));
Fext_d = zeros(Nnodes, 3, length(temps));
Fext = Fext_g+Fext_i+Fext_d;
acc = zeros(length(temps),1);

for s = 1:length(temps);
    t = temps(s);
    uc = vpa(subs(u));
    D = 0.5*Rhoa*uc^2*Cd*S;
    a = g - D/(m_t);
    acc(s) = a;
    
    %Gravity force
for i=1:Nelements
    Fext(Tn(i,1),2,s) = 3; %Actuen tots en la direccio z
    Fext(Tn(i,2),2,s) = 3; 
    Fext(Tn(i,1),1,s) = Tn(i,1);
    Fext(Tn(i,2),1,s) = Tn(i,2);
    Fext(Tn(i,1),3,s) = Fext(Tn(i,1),3,s)+M(i)*g/2;
    Fext(Tn(i,1),3,s) = Fext(Tn(i,1),3,s)-M(i)*a/2;
    Fext(Tn(i,2),3,s) = Fext(Tn(i,2),3,s)+M(i)*g/2;
    Fext(Tn(i,2),3,s) = Fext(Tn(i,2),3,s)-M(i)*a/2;
end
F1 = -D+m_r*(g-a);
Fext(6,3,s) = Fext(6,3,s)+ F1/16;
Fext(8,3,s) = Fext(8,3,s)+ F1/16;
Fext(12,3,s) = Fext(12,3,s)+ F1/16;
Fext(14,3,s) = Fext(14,3,s)+ F1/16;

Fext(7,3,s) = Fext(7,3,s)+ F1/8;
Fext(9,3,s) = Fext(9,3,s)+ F1/8;
Fext(11,3,s) = Fext(11,3,s)+ F1/8;
Fext(13,3,s) = Fext(13,3,s)+ F1/8;

Fext(10,3,s) = Fext(10,3,s)+F1/4;
Fext(1,3,s) = Fext(1,3,s)+m_p*(g-a);
end

    
    
    
    
    
%     
%     
%     
%     Fext_g(Tn(i,1),3,s) = Fext_g(Tn(i,1),3,s)+M(i)*g/2;
%     Fext_g(Tn(i,2),3,s) = Fext_g(Tn(i,2),3,s)+M(i)*g/2;
%     Fext_g(Tn(i,1),1,s) = Tn(i,1);
%     Fext_g(Tn(i,2),1,s) = Tn(i,2);
% end
% 
% % Inertial Forces
% for i=1:Nelements
%     Fext_i(Tn(i,1),3,s) = Fext_i(Tn(i,1),3,s)- M(i)*a/2;
%     Fext_i(Tn(i,2),3,s) = Fext_i(Tn(i,2),3,s)- M(i)*a/2;
% end
% 
% % Drag Forces
% F1 = -D+m_r*g-m_r*a;
% Fext_d(1,3,s) = m_p*(g-a);
% Fext_d(6,3,s) = F1/16;
% Fext_d(7,3,s) = F1/8;
% Fext_d(8,3,s) = F1/16;
% Fext_d(9,3,s) = F1/8;
% Fext_d(10,3,s) = F1/4;
% Fext_d(11,3,s) = F1/8;
% Fext_d(12,3,s) = F1/16;
% Fext_d(13,3,s) = F1/8;
% Fext_d(14,3,s) = F1/16;
% acc(s,1) = a;
% end
% Fext = Fext_g + Fext_i + Fext_d;
% sum(Fext)
% end
