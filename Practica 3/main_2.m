

l_var = 22170;
Nelements = 6*100; NdofsXnode = 2; NnodesXelement = 2; 
L = 36; L1 = 4; L2 = 12; M = 70000; g = 9.81; Me = 2000; E = 70e9;
a = 0.01 ; b = 0.1 ; h = 0.5 ; t = 0.005 ;
Iy = 4.583e-7;
Iz = 1.635e-4;
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
    Kel(:,:,i) = E*Iz/l^3*[12 6*l -12 6*l
        6*l 4*l -6*l 2*l^2
        -12 -6*l 12 -6*l
        6*l 2*l^2 -6*l 4*l^2];
    Nnodes = length(x);
    Ndofs = NdofsXnode*Nnodes;  
    Tn(i,1) = i;
    Tn(i,2) = i+1;
end
%Tenim dos graus de llibertat, el moviment vertical i la rotació 
Td = conncetDOFs(Nelements,NnodesXelement,NdofsXnode,Tn); 
KG = assemblyKG(Nelements,NdofsXnode*NnodesXelement,Ndofs,Td,Kel);
%Reorganitzem la força que tenim per a cada node aquesta matriu contindrà 
%la força i el moment en cada node.
f = zeros(Nnodes, 2); 
for i=1:Nnodes
    if i>=2 && i<=600
    f(i,1) = F(i,1)+F(i-1,3);
    f(i,2) = F(i,2)+F(i-1,4);
    else if i == 1
    f(i,1) = F(i,1);
    f(i,2) = F(i,2);
    else if i == Nnodes
    f(i,1) = F(i-1,3);
    f(i,2) = F(i-1,4);
        end
        end
    end
end
%sumem els moments deguts al motor
for i =1:Nnodes
    f(i,2) = f(i,2)-M*g*(x(i)-L2/2);
end

%Sumem les forces degudes als motor
for i=1:Nelements
    if x(i) == L2/2
        f(i,1) = f(i,1)-Me*g;
    end
end

%Ara tenim la matriu KG, la matriu Td, i la matriu f (Fext).
f_dof = zeros(Nnodes*2,1);
for i=1:Nnodes
    f_dof(i) = f(i,1);
    f_dof(i+1) = f(i,2);
end
% Ara només necessitem la matriu de nodes que volem fixar 
 fixNod =[ 1 1 0
            1 2 0];
        
 [u,R] = solveSys(NdofsXnode,Ndofs,fixNod,KG,f_dof);

