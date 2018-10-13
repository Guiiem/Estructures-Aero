clear all;

%%1. Data input

l_var = 5000; %Lift distribution parameter
Nelements = 6*100; %Number of elements (should be multiple of 6)
NdofsXnode = 2; NnodesXelement = 2; 
L = 36; L1 = 4; L2 = 12; %Geometric data
M = 70000; %Mass of the airplane
g = 9.81; %Gravity acceleration
Me = 2000; %Mass of the motors
E = 70e9; %Young modulus
a = 0.01 ; b = 0.1 ; h = 0.5 ; t = 0.005 ; %Section data
Iy = 4.583e-7; %Inertia y
Iz = 1.635e-4; %Intertia z

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
ver = zeros(Nelements+1); %Vertical displacement
rot = zeros(Nelements+1); %Rotation

for i = 1:Nelements
    %Aixo que hi ha a continuacio que es???
    x_av = (x(i)+x(i+1))/2; %Position of each node?
    l(i) = (x(i+1)-x(i))/2; %Length of each node?
    
    %Lift distribution calculation
    %Te sentit la formula que ens donen????? 
    if x_av<= L1/2
        ql = l_var*(0.85-0.15*cos(2*pi*x(i)/L1));
        q_lift(i) = ql;
    else
        ql = l_var*(L^2-x(i)^2)/(L^2-L1^2);
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
    vector_1(i,:)=[0.5*q(i)*l(i); 1/12*q(i)*l(i)^2; 0.5*q(i)*l(i); -1/12*q(i)*l(i)^2];
    vector_2(i,:)=[1; 1/6*l(i); 1; -1/6*l(i)];
    F_element(i,:) = vector_1(i,:).*0.5*q(i)*l(i).*vector_2(i,:); 
   
    %Punctual force of the motor
    if x(i+1) == 6 %En aquest node ens trobarem la for�a del motor a la dreta
       F_element(i,:) = F_element(i,:)+[0 0 -1/2*Me*g 0];
    end
    if x(i) ==6 %En aquest node ens trobarem la for�a del motor a l'esquerra
       F_element(i,:) = F_element(i,:)+[-1/2*Me*g 0 0 0];
    end
    
    %Calculation of the stifness matrix for each node 
    l = l(i);
    Kel(:,:,i) = E*Iz/l^3*[12 6*l -12 6*l
        6*l 4*l -6*l 2*l^2
        -12 -6*l 12 -6*l
        6*l 2*l^2 -6*l 4*l^2];
    
    %Computation of needed process data
    Nnodes = length(x);
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
    if i>=2 && i<=600
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

%Computation of the moments due to the motor weight

for i =1:m %After the motor location it won't generate any moment
    f_node(i,2) = f_node(i,2)-M*g*(x(i)-L2/2);
end

%Computation of the motor weight force
%No ho haviem afegit ja abans????????????????????
for i=1:Nelements
    if x(i) == L2/2
        f_node(i,1) = f_node(i,1)-Me*g;
    end
end

%Ara tenim la matriu KG, la matriu Td, i la matriu f (Fext).

f_dof = zeros(Nnodes*2,1); %Force on each DOF
for i=1:Nnodes
    f_dof(i) = f_node(i,1);
    f_dof(i+1) = f_node(i,2);
end

% Fixed DOF: The fist node (the center of the fuselage) is completely fixed. 
 fixNod =[ 1 1 0
            1 2 0];
 
 %Solving the system
 [u,R] = solveSys(NdofsXnode,Ndofs,fixNod,KG,f_dof);

 for i=1:Nelements+1
     ver(i) = u(2*i-1);
     rot(i) = u(2*i);
 end
 
 figure;
 plot(x,ver);
 title('despla�ament vertical');
 
 figure;
 plot(x,rot);
 title('rotacio');
 
 

