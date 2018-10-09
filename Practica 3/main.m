%Discretitzarem l'ala en diferents parts. Es realitzarà un mètodes
%semblants a nodes centrats. N nodes i N+1 divisions.
%Dividirem l'ala en la meitat i multipliquem per 2.
%Dades
clear all;
clc;
L = 36; L1 = 4; L2 = 12; M = 70000; g = 9.81; Me = 2000; 

l_var = 22170;
Nnodes = 1000*6;
div = Nnodes+1;
x = linspace(0,L/2, div);
f_lift = zeros(div,1);
F_Lift_node = zeros(Nnodes,1);
f_m = zeros(div,1);
F_motor_node = zeros(Nnodes,1);
f_dis = zeros(div,1);
F_dis_node = zeros(Nnodes, 1);
F_tot = zeros(Nnodes,1);

%% Calcul de les forces
for i=1:div
    if x(i) <= L1/2
        ql = l_var*(0.85-0.15*cos(2*pi*x(i)/L1));
        f_lift(i) = ql;
    else
        ql = l_var*(L^2-x(i)^2)/(L^2-L1^2);
        f_lift(i) = ql;
    end
end
for i=1:Nnodes+1
    if i == Nnodes+1
    F_Lift_node(i) = F_Lift_node(i-1);
    else
    F_Lift_node(i) = (x(i+1)-x(i))*(f_lift(i+1)+f_lift(i))/2; %area trapezi
    end
        
end

%Pel càlcul de la massa dels motors
for i=1:div
if x(i) == 6
    f_m(i) = -M*g;
else
    f_m(i) = 0;
end
end
for i=1:Nnodes+1
    if i == Nnodes+1
    F_motor_node(i) = F_motor_node(i-1);
    else
    F_motor_node(i) = (x(i+1)-x(i))*(f_m(i+1)+f_m(i))/2; %area trapezi
    end
        
end

%Càlcul de la força de la distribució de massa

for i=1:div
    if x(i) <= L1/2
        f_dis(i) = 27/5*M/L*g;
    else
        f_dis(i) = 9/20*M/L*g;
    end
end
for i=1:Nnodes+1
    if i == Nnodes+1
    F_dis_node(i) = F_dis_node(i-1);
    else
        F_dis_node(i) = -(x(i+1)-x(i))*(f_dis(i+1)+f_dis(i))/2;
    end
end
for i=1:Nnodes+1
F_tot(i) = F_dis_node(i)+F_motor_node(i)+F_Lift_node(i);
end


%% Plots de les forces
figure;
plot(x,F_Lift_node); hold on;
plot(-x,F_Lift_node);
title('Força deguda al lift');
figure; 
plot(x,F_dis_node); hold on;
plot(-x,F_dis_node);
title('Distribucio de pes');
figure;
plot(x,F_motor_node); hold on;
plot(-x,F_motor_node);
title('Força dels motors');
figure;
plot(x,F_tot); hold on;
plot(-x,F_tot);
title('Carrega total');
massa_dist = sum(F_dis_node);
massa_motors = sum(F_motor_node);
massa_total = massa_dist + massa_motors;
lift_total = sum(F_Lift_node);
resta= lift_total-massa_total;


%% APARTAT DE Ty

Ty = zeros(Nnodes,1);
for i = Nnodes+1:-1:1
    if i >=2
    Ty(Nnodes+1) = 0;
    Ty(Nnodes) = 0;
    Ty(i-1) = Ty(i)+(F_tot(i)+F_tot(i-1))/((x(i)-x(i-1))*2);
    else
    end
end
figure;
plot(x,Ty); hold on;
plot(-x,Ty)
title('Shear Stress');

Mz = zeros(div, 1);

for i=Nnodes+1:-1:1
   Mz(Nnodes+1) = 0;
   Mz(Nnodes)=0;
   if i>=2
   Mz(i-1) = Mz(i)-(Ty(i)+Ty(i-1))*(x(i)-x(i-1))/2;
   else
   end
end
figure;
plot(x,Mz); hold on;
plot(-x,Mz);
title('bending moment');