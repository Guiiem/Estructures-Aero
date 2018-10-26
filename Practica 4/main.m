clear;
clc;
close all;

%% CALCUL D'INERCIA

S_A = 6.283e-3; %Calculada a mà
Iz_A = 3.1494e-5; %Calculada a mà
b = 60e-3;
t = 80e-3;
h = 8e-3;
Iz_B = 2*((1/12)*b*t^3+b*t*(h/2)^2)+t*(h-t)^3*(1/12);
S_B = (h-t)*t+2*b*t;

%% INPUT DATA
Longitud = 0.7;
Vel = 240/3.6;
mu = 0.4;
H_1 = 1.8;
H_2 = 2.5;
r_roda = 0.22;
alpha = deg2rad(8);
g = 9.81;
E = 200e9;
I_0 = 150; %inercia de la roda

%Definiço de la mtraiu Tn (3 elements)

Tn = [
      1 2
      2 3
      2 4
    ];

%Definició de la matriu de posicions dels diferents nodes (x,y,z)
x = [   
     0    0 0
     0    H_2-H_1 0
     0    H_2 0
     Longitud    H_2 0
    ];

%Definició de la matriu 'mat' amb les propietats del dos materials i
%'Tmat'.

mat = [Iz_A S_A
    Iz_B S_B];

Tmat = [
        1
        1
        2  
    ];

%% Definició de les forces

N = (I_0*Vel)/(r_roda^2*1.5*mu); %força normal de la roda

Fext = [
    1, 1, N*sin(alpha); %Normal no alineada amb la horitzontal
    1, 2, N*cos(alpha);
    1, 3, 0;
    ];

%% PREPROCESS

% Fixem els nodes 4, 3 i 1, de manera que fixem 5 DOF.
fixNod = [% Node        DOF  Magnitude
               4,         1,         0; %Node 4 en la direcció x
               4,         2,         0; %Node 4 en la direcció y
               3,         1,         0; %Node 3 en la direcció x
               3,         2,         0; %Node 3 en la direcció y
               1,         3,         0; %Node 1 en la direcció z (rotació)
];

Ndim = 3;
Nelements = size(Tn, 1);                    % Elements totals 
Nnodes = Nelements+1;                       % Nodes totals (elements +1)
NdofsXnode = 3;                             % Number of DOFs for each node (2 mov + rot)
Ndofs = NdofsXnode*Nnodes;                  % DOFs totals
NnodesXelement = 2;                         % Nodes per cada element
NdofsXelement = NdofsXnode*NnodesXelement;  % Graus de llibertat per cada element

%% RESOLUCIÓ DEL SISTEMA AMB FUNCIONS


Td = connectDOFs(Nelements, NnodesXelement, NdofsXnode, Tn);
Kel = computeKelBar(Nelements, NdofsXelement, x, Tn, Tmat, mat, E);
KG = assemblyKG(Ndofs, NdofsXelement, Nelements, Td, Kel);
f = computeF(NdofsXnode, Ndofs, Fext);
[u,R] = solveSys(NdofsXnode,Ndofs,fixNod,KG,f);
plotDisp(Ndim,Nnodes,u,x,Tn,1);