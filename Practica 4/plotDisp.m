function plotDisp(Ndim,Nnodes,u,x,Tn,fact)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  Ndim     Problem's dimensions
%                  Nnodes   Total number of nodes
%   - u     Global displacement vector [Ndofs x 1]
%            u(I) - Displacement corresponding to the global DOF I
%   - x     Nodal coordinates matrix [Nnodes x Ndim]
%            x(a,i) - Coordinates of node a in the i dimension
%   - Tn    Nodal connectivities table [Nelements x NnodesXelement]
%            Tn(e,a) - Nodal number associated to node a of element e
%   - fact  Scale factor for the deformed structure
%--------------------------------------------------------------------------

% Reshape matrices for plot
U = reshape(u,Ndim,Nnodes);
for i = 1:Ndim
    X0{i} = reshape(x(Tn,i),size(Tn))';
    X{i} = X0{i}+fact*reshape(U(i,Tn),size(Tn))';
end
D = reshape(sqrt(sum(U(:,Tn).^2,1)),size(Tn))';

% Open and initialize figure
figure('color','w');
hold on;       % Allow multiple plots on same axes
box on;        % Closed box axes
axis equal;    % Keep aspect ratio to 1
colormap jet;  % Set colormap colors
set(axes, 'Xdir', 'reverse');
set(gca,'xtick',[],'ytick',[],'units','pixels'); % Take out ticks from axes

% Plot undeformed structure
patch(X0{:},zeros(size(D)),'edgecolor',[0.5,0.5,0.5],'linewidth',2);

% Plot deformed structure with displacement magnitude coloring
patch(X{:},D,'edgecolor','interp','linewidth',2);

% Set colorbar properties
caxis([min(D(:)),max(D(:))]); % Colorbar limits
cbar = colorbar;              % Create colorbar
set(cbar,'Ticks',linspace(min(D(:)),max(D(:)),5))

% Set title
title('Displacement');