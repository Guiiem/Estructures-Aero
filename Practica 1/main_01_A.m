%-------------------------------------------------------------------------%
% ASSIGNMENT 01-A
%-------------------------------------------------------------------------%
% Date:
% Author/s:
%

clear;
close all;

%% INPUT DATA

F = 100;
E = 5e10;
A = 80e-6;

%% PREPROCESS

% Nodal coordinates matrix creation
%  x(a,j) = coordinate of node a in the dimension j
x = [%     X      Y
           0,     0; % (1)
         500,   200; % (2)
        1000,   400; % (3)
        1500,   600; % (4)
           0,   500; % (5)
         500,   600; % (6)
        1000,   700; % (7)
        1500,   800; % (8)
]/1e3;

% Connectivities matrix ceation
%  Tn(e,a) = global nodal number associated to node a of element e
Tn = [%     a      b
           1,     2; % (1)
           5,     2; % (2)
           5,     6; % (3)
           2,     6; % (4)
           2,     3; % (5)
           2,     7; % (6)
           6,     7; % (7)
           3,     7; % (8)
           3,     4; % (9)
           7,     4; % (10)
           7,     8; % (11)
           4,     8; % (12)
];

% External force matrix creation
%  Fext(k,1) = node at which the force is applied
%  Fext(k,2) = DOF (direction) at which the force is applied
%  Fext(k,3) = force magnitude in the corresponding DOF
Fext = [%   Node        DOF  Magnitude   
               2,         2,         F;
               3,         2,         F;
               4,         2,         F;
];

% Fix nodes matrix creation
%  fixNod(k,1) = node at which some DOF is prescribed
%  fixNod(k,2) = DOF prescribed
%  fixNod(k,3) = prescribed displacement in the corresponding DOF (0 for fixed)
fixNod = [% Node        DOF  Magnitude
               1,         1,         0;
               1,         2,         0;
               5,         1,         0;
               5,         2,         0;
];

% Material data
%  mat(m,1) = Young modulus of material m
%  mat(m,2) = Section area of material m
%  --more columns can be added for additional material properties--
mat = [% Young M.   Section A.     
                E,           A;  % Material (1)
];

% Material connectivities
%  Tmat(e) = Row in mat corresponding to the material associated to element e 
Tmat = [% Mat. index
                   1; % (1)
                   1; % (2)
                   1; % (3)
                   1; % (4)
                   1; % (5)
                   1; % (6)
                   1; % (7)
                   1; % (8)
                   1; % (9)
                   1; % (10)
                   1; % (11)
                   1; % (12)
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

% Computation of element stiffness matrices
Kel = computeKelBar(Ndim,Nelements,x,Tn,mat,Tmat);

% Global matrix assembly
KG = assemblyKG(Nelements,NdofsXelement,Ndofs,Td,Kel);

% Global force vector assembly
f = computeF(NdofsXnode,Ndofs,Fext);

% System resolution
[u,R] = solveSys(NdofsXnode,Ndofs,fixNod,KG,f);

% Compute strain and stresses
[eps,sig] = computeStrainStressBar(Ndim,Nelements,u,Td,x,Tn,mat,Tmat);

%% POSTPROCESS

% Plot displacements
plotDisp(Ndim,Nnodes,u,x,Tn,1);

% Plot strains
plotStrainStress(Ndim,eps,x,Tn);

% Plot stress
plotStrainStress(Ndim,sig,x,Tn);