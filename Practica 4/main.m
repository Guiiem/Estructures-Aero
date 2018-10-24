clear;
clc;
close all;

%% INERTIA
[Iz_B, S_B] = computeInertia_B(60e-3, 80e-3, 8e-3);
Iz_A = 3.1494e-5;
S_A = 6.283e-3;

%% INPUT DATA
Longitud = 0.7;
H1 = 1.8;
H2 = 2.5;
rw = 0.22;
alpha = deg2rad(8);
g = 9.81;
E = 200e9;
I_0 = 150;
V = 240/3.6;
mu = 0.4;

x = [%x   %y
     0    0 0%1
     0    H2-H1 0%2
     0    H2 0%3
     Longitud    H2 0%4
    ];
Tn = [
      1 2;
      2 3;
      2 4;
    ];

mat = [Iz_A S_A
    Iz_B S_B];
Tmat = [
        1
        1
        2  
    ];
N = (I_0*V)/(rw^2*1.5*mu);

Fext = [
    1, 1, N/sin(alpha); %Normal no alineada amb la horitzontal
    1, 2, N/cos(alpha);
    1, 3, 0;
    ];

%% PREPROCESS

% Fixem els nodes 4, 3 i 1, de manera que fixem 5 DOF.
fixNod = [% Node        DOF  Magnitude
               4,         1,         0;
               4,         2,         0;
               3,         1,         0;
               3,         2,         0;
               1,         3,         0;
];

NdofsXnode = 3;                             % Number of DOFs for each node
Ndim = 3;
Nelements = size(Tn, 1);                    % Total number of elements (it needs to be multiple of 6)
Nnodes = Nelements+1;                       % Total number of nodes
Ndofs = NdofsXnode*Nnodes;                  % Total number of degrees of freedom
NnodesXelement = 2;                         % Number of nodes for each element
NdofsXelement = NdofsXnode*NnodesXelement;  % Number of DOFs for each element

%% RESOLUCIÓ DEL SISTEMA AMB FUNCIONS


Td = connectDOFs(Nelements, NnodesXelement, NdofsXnode, Tn);
Kel = computeKelBar(Nelements, NdofsXelement, x, Tn, Tmat, mat, E);
KG = assemblyKG(Ndofs, NdofsXelement, Nelements, Td, Kel);
f = computeF(NdofsXnode, Ndofs, Fext);
[u,R] = solveSys(NdofsXnode,Ndofs,fixNod,KG,f);
plotDisp(Ndim,Nnodes,u,x,Tn,1);