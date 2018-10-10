

l_var = 22170;
Nelements = 6*100;
L = 36; L1 = 4; L2 = 12; M = 70000; g = 9.81; Me = 2000; E = 70e9;
a = 0.01 ; b = 0.1 ; h = 0.5 ; t = 0.005 ;
Iy = inercia(b,h,t,a);
Iz = 1/12*2*t*b;
x = linspace(0,L/2, Nelements+1);
q_lift = zeros(Nelements,1);
q_dis = zeros(Nelements,1);
q = zeros(Nelements,1);
l = zeros(Nelements,1);
F = zeros(Nelements,4);
vector_1 = zeros(Nelements,4);
vector_2 = zeros(Nelements,4);
Kel = zeros(4, 4, Nelements);
for i = 1:Nelements
    x_av = (x(i)+x(i+1))/2;
    l(i) = (x(i+1)-x(i))/2;
    %Calcul de la força del lift
    if x_av<= L1/2
        ql = l_var*(0.85-0.15*cos(2*pi*x(i)/L1));
        q_lift(i) = ql;
    else
        ql = l_var*(L^2-x(i)^2)/(L^2-L1^2);
        q_lift(i) = ql;
    end
    
    %calcul de la força del pes
    if x_av <= L1/2
        q_dis(i) = 27/5*M/L*g;
    else
        q_dis(i) = 9/20*M/L*g;
    end
    q(i) = q_dis(i)+q_lift(i); %Sumem les dos forces per tenir totes les distribucions
    vector_1(i,:)=[0.5*q(i)*l(i); 1/12*q(i)*l(i)^2; 0.5*q(i)*l(i); -1/12*q(i)*l(i)^2];
    vector_2(i,:)=[1; 1/6*l(i); 1; -1/6*l(i)];
    F(i,:) = vector_1(i,:).*0.5*q(i)*l(i).*vector_2(i,:);
   
    if x(i+1) == 6 %En aquest node ens trobarem la força del motor a la dreta
       F(i,:) = F(i,:)+[0 0 -1/2*Me*g 0];
    end
    if x(i) ==6 %En aquest node ens trobarem la força del motor a l'esquerra
       F(i,:) = F(i,:)+[-1/2*Me*g 0 0 0];
    end
    l = l(i);
    Kel(:,:,i) = [12 6^l -12 6*l
        6*l 4*l -6*l 2*l^2
        -12 -6*l 12 -6*l
        6*l 2*l^2 -6*l 4*l^2];
   
    
     %Ara calculem l'equació de cada element
end

    

