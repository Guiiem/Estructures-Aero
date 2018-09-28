function Kel = computeKelBar(Ndim,Nelements,x,Tn,mat,Tmat)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  Ndim        Problem's dimensions
%                  Nelements   Total number of elements
%   - x     Nodal coordinates matrix [Nnodes x Ndim]
%            x(a,i) - Coordinates of node a in the i dimension
%   - Tn    Nodal connectivities table [Nelements x NnodesXelement]
%            Tn(e,a) - Nodal number associated to node a of element e
%   - mat   Material properties table [Nmat x NpropertiesXmat]
%            mat(m,1) - Young modulus of material m
%            mat(m,2) - Section area of material m
%   - Tmat  Material connectivities table [Nelements]
%            Tmat(e) - Material index of element e
%--------------------------------------------------------------------------
% It must provide as output:
%   - Kel   Elemental stiffness matrices [NdofsXelement x NdofsXelement x Nelements]
%            Kel(i,j,e) - Term in (i,j) position of stiffness matrix for element e
%--------------------------------------------------------------------------

% Longituds de barra
j=0;
for j=1:Nelements
    if(Ndim==2) % 2D
   L(j,1)=[sqrt((x(Tn(j,1),1)-(x(Tn(j,2),1)))^2 + ...
       + (x(Tn(j,1),2)-(x(Tn(j,2),2)))^2)];
    else % 3D
    L(j,1)=[sqrt((x(Tn(j,1),1)-(x(Tn(j,2),1))^2 + ...
       + (x(Tn(j,1),2)-(x(Tn(j,2),2)))^2 + ...
       + (x(Tn(j,1),3)-(x(Tn(j,2),3)))^2))];
    end
end    

% Matriu R per cada element
j=0; 
    for j=1:Nelements
        if(Ndim==2) % 2D
            % fila x col x profunditat 
R(1,1,j)= 1/L(j)*( (x(Tn(j,2),1)-x(Tn(j,1),1))); R(2,3,j)=1/L(j)*( (x(Tn(j,2),1)-x(Tn(j,1),1)));    
R(1,2,j)= 1/L(j)*(x(Tn(j,2),2)-x(Tn(j,1),2)); R(2,4,j)= 1/L(j)*(x(Tn(j,2),2)-x(Tn(j,1),2));
R(1,3)=0; R(1,4)=0; R(2,1)=0; R(2,2)=0;

        else % 3D
R(1,1,j)= 1/L(j)*( (x(Tn(j,2),1)-x(Tn(j,1),1))); 
R(1,2,j)= 1/L(j)*(x(Tn(j,2),2)-x(Tn(j,1),2));
R(1,3,j)= 1/L(j)*(x(Tn(j,2),3)-x(Tn(j,1),3));

R(2,4,j)=1/L(j)*( (x(Tn(j,2),1)-x(Tn(j,1),1)));    
R(2,5,j)= 1/L(j)*(x(Tn(j,2),2)-x(Tn(j,1),2));
R(2,6,j)= 1/L(j)*(x(Tn(j,2),2)-x(Tn(j,1),2)); 

R(1,4)=0; R(1,5)=0; R(1,6)=0; R(2,1)=0; R(2,2)=0; R(2,3)=0;

        end
    end
    
% K prima
e=0;
 for e=1:Nelements
        % fila x col x profunditat
   k=mat(Tmat(e),1)*mat(Tmat(e),2)/L(e,1)*1; %E*A/L
        %if(Ndim==2) % 2D
   Kprima(1,1,e) = k; Kprima(1,2,e) = -k;
   Kprima(2,1,e) = -k; Kprima(2,2,e) = k; %Sempre es la mateixa matriu 2x2 WHYYY?
        %else % 3D
   %Kprima(1,1,e) = k; Kprima(1,2,e) = -k; Kprima(1,3,e) = -k;
   %Kprima(2,1,e) = -k; Kprima(2,2,e) = k; Kprima(2,3,e) = -k;    
   %Kprima(3,1,e) = -k; Kprima(3,2,e) = -k; Kprima(3,3,e) = k;         
        %end
 end
 
j=0;
for j=1:Nelements
    Kel(:,:,j) = R(:,:,j).'*Kprima(:,:,j)*R(:,:,j);
end

end
