%-------------------------------------------------------------------------%
% ASSIGNMENT 01-A
%-------------------------------------------------------------------------%
% Date: 
% Author/s: Guillem Verg�s Plaza; David Rodriguez Pozo.
%

clear all;
close all;
clc;
%% INPUT DATA

Db = 7.8e-3;
db = 3.8e-3;
Dc = 1.5e-3;
Ec = 200000e6; %Young's modulus-bars
Eb = 70000e6; %Young's modulus-cables
Ac = 1.77e-6; %Cable area
Ab = 36.44e-6; %Bars area
Rhoa = 1.225; %Air density
Rhoc = 1500;  %Cable density
Rhob = 2300; %Bars density
Ib = 1/4*pi*((Db/2)^4-(db/2)^4);
Ic = 1/4*pi*(Dc/2)^4;
S = 17.5; %Parachute surface
Cd = 1.25; %Drag coeficient
t_s = 0.001; %Thickness of the robe
rhos = 1500; %Desnity of the robe
g = 9.81; %Gravity force
temps = 0:2:10;
m_p = 120; 
m_r = S*t_s*rhos;

%% PREPROCESS

% Nodal coordinates matrix creation
%  x(a,j) = coordinate of node a in the dimension j
x = [%   X       Y       Z
     0.000,  0.000,  0.000; % 1
    -0.725, -0.731,  3.250; % 2
     0.725, -0.731,  3.250; % 3
     0.725,  0.731,  3.250; % 4
    -0.725,  0.731,  3.250; % 5
     2.900, -1.450,  5.630; % 6
     2.900,  0.000,  5.820; % 7
     2.900,  1.450,  5.630; % 8
     0.000,  1.450,  6.340; % 9
     0.000,  0.000,  6.530; % 10
     0.000, -1.450,  6.340; % 11
    -2.900,  1.450,  5.630; % 12
    -2.900,  0.000,  5.820; % 13
    -2.900, -1.450,  5.630; % 14
];

% Connectivities matrix ceation
%  Tn(e,a) = global nodal number associated to node a of element e
Tn = [%   A     B
         1     2  % 1
         1     3  % 2
         1     4  % 3
         1     5  % 4
         2     3  % 5
         3     4  % 6 
         4     5  % 7 
         5     2  % 8
         2     4  % 9
         2    10  % 10
         2    11  % 11
         2    13  % 12
         2    14  % 13
         3     6  % 14
         3     7  % 15
         3    10  % 16
         3    11  % 17
         4     7  % 18
         4     8  % 19
         4     9  % 20
         4    10  % 21
         5     9  % 22
         5    10  % 23
         5    12  % 24
         5    13  % 25
        14    11  % 26
        11    10  % 27
        14    10  % 28
        14    13  % 29
        13    10  % 30
        13    12  % 31
        12    10  % 32
        12     9  % 33
         9    10  % 34
         9     8  % 35
         8     7  % 36
        10     7  % 37
        10     8  % 38
        11     6  % 39
         6     7  % 40
        10     6  % 41
];


% External force matrix creation
%  Fext(k,1) = node at which the force is applied
%  Fext(k,2) = DOF (direction) at which the force is applied
%  Fext(k,3) = force magnitude in the corresponding DOF



% Fix nodes matrix creation
%  fixNod(k,1) = node at which some DOF is prescribed
%  fixNod(k,2) = DOF prescribed
%  fixNod(k,3) = prescribed displacement in the corresponding DOF (0 for fixed)

fixNod = [ % Node DOF Magnitude
    1 , 1 , 0 ; %El punt 1 queda fix i evitem translaci�.
    1 , 2 , 0 ;
    1 , 3 , 0 ;
    2 , 1 , 0 ; %Evitem que el punt 10 es mogui en x i y.
    2 , 2 , 0 ;
    3 , 2 , 0 ; %Evitem que el punt 7 es mogui en la direccio y.
];
% Material data
%  mat(m,1) = Young modulus of material m
%  mat(m,2) = Section area of material m
%  --more columns can be added for additional material properties--
mat = [% Young M.   Section A.     
                Ec, Ac, Rhoc, Ic;  % Cables
                Eb, Ab, Rhob, Ib; %Bars
                
];

% Material connectivities
%  Tmat(e) = Row in mat corresponding to the material associated to element e 
Tmat = [% Mat. index
         1  % 1
         1  % 2
         1  % 3
         1  % 4
         2  % 5
         2  % 6 
         2  % 7 
         2  % 8
         2  % 9
         1  % 10
         1  % 11
         1  % 12
         1  % 13
         1  % 14
         1  % 15
         1  % 16
         1  % 17
         1  % 18
         1  % 19
         1  % 20
         1  % 21
         1  % 22
         1  % 23
         1  % 24
         1  % 25
         2  % 26
         2  % 27
         2  % 28
         2  % 29
         2  % 30
         2  % 31
         2  % 32
         2  % 33
         2  % 34
         2  % 35
         2  % 36
         2  % 37
         2  % 38
         2  % 39
         2  % 40
         2  % 41
];

%% SOLVER

% Dimensions
Ndim = size(x,2);                           % Number of dimensions
NdofsXnode = Ndim;                          % Number of DOFs for each node
Nnodes = size(x,1);                         % Total number of nodes
Ndofs = NdofsXnode*Nnodes;                  % Total number of degrees of freedom
Nelements = size(Tn,1);                     % Total number of elements
NnodesXelement = size(Tn,2);                % Number of nodes for each element
NdofsXelement = NdofsXnode*NnodesXelement;  % Number of DOFs for each element 

% Computation of the DOFs connectivities
Td = connectDOFs(Nelements,NnodesXelement,NdofsXnode,Tn);

% Computation of element stiffness matrices and the length and mass 
[Kel, M] = computeKelBar(Ndim,Nelements,x,Tn,mat,Tmat);

% Global matrix assembly
KG = assemblyKG(Nelements,NdofsXelement,Ndofs,Td,Kel);
massa_bc = sum(M);
m_total = m_r+massa_bc+m_p;

% Compute the velocity as a function of time
u = solveEDO (m_total, Cd, Rhoa,S, g);


% Fext matrix computation
[Fext,acc] = Fext(Nnodes, M, u, temps, Rhoa, Cd, g, m_p, m_r, Nelements, Tn, m_total, S);

% Global force vector assembly
f = computeF(NdofsXnode,Ndofs,Fext);

% System resolution
[U,R] = solveSys(NdofsXnode,Ndofs,fixNod,KG,f);

% Compute strain and stresses
[eps,sig,eps_t,sig_t] = computeStrainStressBar(Ndim,Nelements,U,Td,x,Tn,mat,Tmat,temps);

Udef = U(:,:,length(temps));

%Conseguim els maxims i minims
[maxsig, minsig] = getmax(sig_t);

%Buckling 
SigLim = zeros(Nelements,1);
SigLim(:,1)=sig_t(:,1,length(temps)); %Agafo dels ultims valors pq tenen els valors de compressio mes grans
sigcr = sigvinc(Nelements, x, Ic, Ib, mat, Tmat, Tn);
%Vinclament: Compressio
   for j=1:41
       if (SigLim(j,1)<0) %Nomes ens interessen els elements a compressio
       if (abs(SigLim(j,1))>=abs(sigcr(j,1)))
           BUCK(j,1)=1;
       else
           BUCK(j,1)=0;
       end
       else 
           BUCK(j,1)=0;
       end
   end




%% POSTPROCESS
%Plot max and min. 
plot(temps,maxsig);
figure();
plot(temps,minsig);

% Plot displacements
plotDisp(Ndim,Nnodes,Udef,x,Tn,1);

% Plot strains
plotStrainStress(Ndim,eps,x,Tn);

% Plot stress
plotStrainStress(Ndim,sig,x,Tn);

%Plot pandeo
plotStrainStress(Ndim,BUCK, x, Tn);