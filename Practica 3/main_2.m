

l_var = 22170;
Nelements = 6*100;
L = 36; L1 = 4; L2 = 12; M = 70000; g = 9.81; Me = 2000;
x = linspace(0,L/2, Nelements+1);
q_lift = zeros(Nelements,1);
q_dis = zeros(Nelements,1);
q = zeros(Nelements,1);
l = zeros(Nelements,1);
F = zeros(Nelements,1);
vector_1 = zeros(Nelements,4);
vector_2 = zeros(Nelements,4);
for i = 1:Nelements
    x_av = (x(i)+x(i+1))/2;
    l(i) = (x(i+1)-x(i))/2;
    %Calcul de la for�a del lift
    if x_av<= L1/2
        ql = l_var*(0.85-0.15*cos(2*pi*x(i)/L1));
        q_lift(i) = ql;
    else
        ql = l_var*(L^2-x(i)^2)/(L^2-L1^2);
        q_lift(i) = ql;
    end
    
    %calcul de la for�a del pes
    if x_av <= L1/2
        q_dis(i) = 27/5*M/L*g;
    else
        q_dis(i) = 9/20*M/L*g;
    end
    q(i) = q_dis(i)+q_lift(i); %Sumem les dos forces per tenir totes les distribucions
    vector_1(i,:)=[0.5*q(i)*l(i); 1/12*q(i)*l(i)^2; 0.5*q(i)*l(i); -1/12*q(i)*l(i)^2];
    vector_2(i,:)=[1; 1/6*l(i); 1; -1/6*l(i)];
    F(i) = vector_1(i,:).*0.5*q(i)*l(i).*vector_2(i,:);
    %Ara calculem l'equaci� de cada element
    
     
end

    

