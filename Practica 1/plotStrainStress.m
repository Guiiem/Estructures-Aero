function plotStrainStress(Ndim,s,x,Tn)
%--------------------------------------------------------------------------
% The function takes as inputs:
%   - Dimensions:  Ndim            Problem's dimensions
%   - s     Strain/stress vector [Nelements x 1]
%            s(e) - Strain/stress associated to bar e
%   - x     Nodal coordinates matrix [Nnodes x Ndim]
%            x(a,i) - Coordinates of node a in the i dimension
%   - Tn    Nodal connectivities table [Nelements x NnodesXelement]
%            Tn(e,a) - Nodal number associated to node a of element e
%   - fact  Scale factor for the deformed structure
%--------------------------------------------------------------------------

% Reshape matrices for plot
for i = 1:Ndim
    X0{i} = reshape(x(Tn,i),size(Tn))';
end
S = repmat(s',size(Tn,2),1);

% Open and initialize figure
figure('color','w');
hold on;       % Allow multiple plots on same axes
box on;        % Closed box axes
axis equal;    % Keep aspect ratio to 1
colormap jet;  % Set colormap colors
set(gca,'xtick',[],'ytick',[],'units','pixels'); % Take out ticks from axes

% Plot undeformed structure
patch(X0{:},S,'edgecolor','flat','linewidth',2);

% Set colorbar properties
caxis([min(S(:)),max(S(:))]);
cbar = colorbar;
set(cbar,'Ticks',linspace(min(S(:)),max(S(:)),5));