

l_var = 22170;
Nelements = 6*100;
x = linspace(0,L/2, Nelements+1);
L = 36; L1 = 4; L2 = 12; M = 70000; g = 9.81; Me = 2000; 
f_lift = zeros(Nelements,1);
f_m = zeros(Nelements,1);
f_dis = zeros(Nelements,1);
for i = 1:Nelements
    x_av = (x(i)+x(i+1))/2;
    
    %Calcul de la força del lift
    if x_av<= L1/2
        ql = l_var*(0.85-0.15*cos(2*pi*x(i)/L1));
        f_lift(i) = ql;
    else
        ql = l_var*(L^2-x(i)^2)/(L^2-L1^2);
        f_lift(i) = ql;
    end
    
    if x(i) <= L1/2
        f_dis(i) = 27/5*M/L*g;
    else
        f_dis(i) = 9/20*M/L*g;
    end
     
end

if x(i) <= L1/2
        f_dis(i) = 27/5*M/L*g;
    else
        f_dis(i) = 9/20*M/L*g;
 end
    

