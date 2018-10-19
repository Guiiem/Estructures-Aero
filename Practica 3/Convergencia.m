clear all;
%%Convergence Study
%estudiarem de 1 a 600 nodes i guardarem els resultats de Ty i de Mz a la
%punta de l'ala

Nvar = 0:6:1800;
Nvar(1) = 1; %perque comenxi per 1;
Mz_max = zeros(length(Nvar));
Ty_max = zeros(length(Nvar));

for j = 1:length(Nvar)
    
Nelements = Nvar(j);

%%1. Data input

l_var = 30455; %Lift distribution parameter %EL BO ES APROX 30455
%Nelements = 6*10; %Number of elements (should be multiple of 6)
NdofsXnode = 2; NnodesXelement = 2; 
L = 36; L1 = 4; L2 = 12; %Geometric data
M = 70000; %Mass of the airplane
g = 9.81; %Gravity acceleration
Me = 2000; %Mass of the motors
E = 70e9; %Young modulus
a = 0.01 ; b = 0.1 ; h = 0.5 ; t = 0.005 ; %Section data
Iy = 4.583e-7; %Inertia y
Iz = 1.635e-4; %Intertia z
Tn = zeros(Nelements,2);

%%2. Matrix creation

x = linspace(0,L/2, Nelements+1); %Discretitation 

q_lift = zeros(Nelements,1); %Aerodynamic force distribution
q_dis = zeros(Nelements,1); %Mass force distribution
q = zeros(Nelements,1); %Total force distribution
l = zeros(Nelements,1);
F_element = zeros(Nelements,4);
vector_1 = zeros(Nelements,4);
vector_2 = zeros(Nelements,4);
Kel = zeros(4, 4, Nelements); %Stifness matrix
ver = zeros(Nelements+1,1); %Vertical displacement
rot = zeros(Nelements+1,1); %Rotation
x_av_vector = zeros(Nelements,1);

for i = 1:Nelements
    
    x_av_vector(i) = (x(i)+x(i+1))/2; %Position of the center of each element
    l(i) = (x(i+1)-x(i)); %Length of each element
    x_av = x_av_vector(i);
    %Lift distribution calculation
    
    if x_av <= L1/2
        ql = l_var*(0.85-0.15*cos(2*pi*x(i)/L1));
        q_lift(i) = ql;
    else
        ql = l_var*((L)^2-4*x(i)^2)/(L^2-L1^2);
        q_lift(i) = ql;
    end
    
    %Mass density distribution calculation
    if x_av <= L1/2
        q_dis(i) = -27/5*M/L*g;
    else
        q_dis(i) = -9/20*M/L*g;
    end
    
    q(i) = q_dis(i)+q_lift(i); %Total force distribution
    
    %Force distribution to ELEMENT discrete forces
    vector_2(i,:)=[1; 1/6*l(i); 1; -1/6*l(i)];
    F_element(i,:) = 0.5*q(i)*l(i)*vector_2(i,:); 
   
    %Punctual force of the motor
    if x(i+1) == 6 %En aquest node ens trobarem la força del motor a la dreta
       F_element(i,:) = F_element(i,:)+[0 0 -1/2*Me*g 0];
    end
    if x(i) == 6 %En aquest node ens trobarem la força del motor a l'esquerra
       F_element(i,:) = F_element(i,:)+[-1/2*Me*g 0 0 0];
    end
    
    %Calculation of the stifness matrix for each node 
    l = l(i);
    Kel(:,:,i) = E*Iz/l^3*[12 6*l -12 6*l
        6*l 4*l^2 -6*l 2*l^2
        -12 -6*l 12 -6*l
        6*l 2*l^2 -6*l 4*l^2];
    %Computation of needed process data
    Nnodes = Nelements+1;
    Ndofs = NdofsXnode*Nnodes;  
    
    %Node connectivities matrix
    Tn(i,1) = i;
    Tn(i,2) = i+1;
    
    %Finding of the motor node
    if x(i) == L2/2
        m = i;  
    end
end


%Two DOF: vertical displacement and rotation of each node

Td = conncetDOFs(Nelements,NnodesXelement,NdofsXnode,Tn); 
KG = assemblyKG(Nelements,NdofsXnode*NnodesXelement,Ndofs,Td,Kel);

%Reorganization of the force for each node: creation of a matrix that
%will have the force and moment of each node.

f_node = zeros(Nnodes, 2); 
for i=1:Nnodes
    if i>=2 && i<=Nelements
    f_node(i,1) = F_element(i,1)+F_element(i-1,3);
    f_node(i,2) = F_element(i,2)+F_element(i-1,4);
    else if i == 1
    f_node(i,1) = F_element(i,1);
    f_node(i,2) = F_element(i,2);
    else if i == Nnodes
    f_node(i,1) = F_element(i-1,3);
    f_node(i,2) = F_element(i-1,4);
        end
        end
    end
end


f_dof = zeros(Nnodes*2,1); %Force on each DOF
for i=1:Nnodes
    f_dof(2*i-1) = f_node(i,1);
    f_dof(2*i) = f_node(i,2);
end

% Fixed DOF: The fist node (the center of the fuselage) is completely fixed. 
 fixNod =[ 1 1 0
            1 2 0];
 
 %Solving the system
 [u,R] = solveSys(NdofsXnode,Ndofs,fixNod,KG,f_dof);

%% Càlcul de esfuerzos cortantes i moments.
Ty = zeros(Nelements,1);
Mz = zeros(Nelements,1);
Ty(1) = -R(1);
Mz(1) = -R(2);
for i = 2:Nelements
    length = abs(x(i-1)-x(i));
    Ty(i) = Ty(i-1)-f_dof(2*i-1);
    Mz(i) = Mz(i-1)-f_dof(2*i)-((Ty(i)+Ty(i-1))/2*length);
end
 
 for i=1:Nelements+1
     ver(i) = u(2*i-1);
     rot(i) = u(2*i);
 end
 
 Mz_max(j) = max(Mz);
 Ty_max(j) = max(Ty);

 
end
 %% Plot de totes les variables necessàries
 
 figure;
 semilogx(Nvar, Mz_max); hold on;
 semilogx(Nvar, Ty_max);
 title('convergencia');
 
 figure;
 plot(x,ver); hold on;
 plot(-x,ver);
 title('desplaçament vertical');
 
 figure;
 plot(x,rot); hold on;
 plot(-x,rot);
 title('rotacio');
 
 figure
 plot(x_av_vector, q_lift); hold on;
 plot(-x_av_vector, q_lift);
 title('lift');
 
 figure;
 plot(x_av_vector, q_dis); hold on;
 plot(-x_av_vector, q_dis);
  plot(x_av_vector, q_lift); hold on;
 plot(-x_av_vector, q_lift);
 title('massa');
 
 figure; 
 plot(x_av_vector, Ty); hold on;
 plot(-x_av_vector, Ty)
 title('Shear force');
 
 figure;
 plot(x_av_vector, Mz); hold on;
 plot(-x_av_vector,Mz);
 title('Bending moment');
