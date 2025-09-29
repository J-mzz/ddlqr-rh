% clear all; clc
%% compare consistency set

% load data and get vertices 
load ddlqr_ic data
data_ic = data;

load ddlqr_rh data
data_rh = data;

clear data

% get vertices from data
vert_ic = getVertices(data_ic);
vert_rh = getVertices(data_rh);

% visualized
% option 1: chosen coordinates
%        2: 2d pca
%        3: 3d pca
opt = 3; 
visualize(vert_ic, vert_rh, opt)

ax = gca;
ax.Title.FontSize = 16;
ax.Legend.FontSize = 16;

view(-31.2302, 13.7764)


%% Local Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function visualize(V1, V2, option)

% Option 1: chosen 3 coordinates
if option == 1
% === Step 1: Take first 3 coordinates directly ===
V1_proj = V1(:,4:6);
V2_proj = V2(:,4:6);

% === Step 2: Compute convex hulls in 3D ===
K1 = convhull(V1_proj(:,1), V1_proj(:,2), V1_proj(:,3));
K2 = convhull(V2_proj(:,1), V2_proj(:,2), V2_proj(:,3));

% === Step 3: Plot ===
figure; hold on; axis equal; grid on

trisurf(K1, V1_proj(:,1), V1_proj(:,2), V1_proj(:,3), ...
        'FaceColor',[0.6 0.7 1], 'FaceAlpha',0.2, 'EdgeColor','b');

trisurf(K2, V2_proj(:,1), V2_proj(:,2), V2_proj(:,3), ...
        'FaceColor',[1 0.6 0.6], 'FaceAlpha',0.5, 'EdgeColor','r');

xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
legend('ours-fixed $\mathcal{C}$','ours-updating $\mathcal{C}$',...
    'Interpreter', 'latex');
title('Consistency sets using first 3 coordinates');

view(3); % default 3D view
camlight; lighting gouraud;

% Option 2: 2D PCA
elseif option == 2
% === Step 1: PCA projection ===
allV = [V1; V2];   % combine vertices
[coeff, ~] = pca(allV);

% Take first 2 principal components
V1_proj = V1 * coeff(:,1:2);
V2_proj = V2 * coeff(:,1:2);

% === Step 2: Convex hulls in 2D ===
K1 = convhull(V1_proj(:,1), V1_proj(:,2));
K2 = convhull(V2_proj(:,1), V2_proj(:,2));

% === Step 3: Compute areas ===
A1 = polyarea(V1_proj(K1,1), V1_proj(K1,2)); % area outer
A2 = polyarea(V2_proj(K2,1), V2_proj(K2,2)); % area inner
shrinkage_ratio = A2 / A1;

% === Step 4: Plot overlay ===
figure; hold on; axis equal; grid on

fill(V1_proj(K1,1), V1_proj(K1,2), [0.8 0.8 1], ...
     'FaceAlpha',0.3, 'EdgeColor','b', 'LineWidth',1.5);

fill(V2_proj(K2,1), V2_proj(K2,2), [1 0.6 0.6], ...
     'FaceAlpha',0.5, 'EdgeColor','r', 'LineWidth',1.5);

xlabel('PC 1'); ylabel('PC 2');
legend('ours-fixed $\mathcal{C}$','ours-updating $\mathcal{C}$',...
    'Interpreter', 'latex');
title(sprintf('Projection of polytopes (Shrinkage ratio = %.3f)', shrinkage_ratio));


% Option 3: 3D PCA
elseif option == 3
allV = [V1; V2];   % combine vertices
[coeff, ~] = pca(allV);

% Project to first 3 principal components
V1_proj = V1 * coeff(:,1:3);
V2_proj = V2 * coeff(:,1:3);

% === Step 2: Convex hulls in 3D ===
K1 = convhull(V1_proj(:,1), V1_proj(:,2), V1_proj(:,3));
K2 = convhull(V2_proj(:,1), V2_proj(:,2), V2_proj(:,3));

% === Step 3: Plot overlay ===
figure; hold on; axis equal; grid on

trisurf(K1, V1_proj(:,1), V1_proj(:,2), V1_proj(:,3), ...
        'FaceColor',[0.6 0.7 1], 'FaceAlpha',0.2, 'EdgeColor','b');

trisurf(K2, V2_proj(:,1), V2_proj(:,2), V2_proj(:,3), ...
        'FaceColor',[1 0.6 0.6], 'FaceAlpha',0.5, 'EdgeColor','r');

xlabel('PC 1'); ylabel('PC 2'); zlabel('PC 3');
legend('ours-fixed $\mathcal{C}$','ours-updating $\mathcal{C}$',...
    'Interpreter', 'latex');
title('Consistency sets projected to first 3 PCA components');

view(3); % set default 3D view
camlight; lighting gouraud; % nicer shading

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vert = getVertices(data)
% get vertices from data = (xk+1, xk, uk)

Xp = data.Xp;
Xm = data.Xm;
U = data.U;
n = size(Xp,2);
eps = data.eps;

aa = []; bb = [];
for i = 1:n
    aa = [aa; kron(eye(2),Xm(:, i)')];
    bb = [bb; kron(eye(2),U(i))];
end
xi = Xp(:);

N = [aa bb; -aa -bb];
if rank(N) ~= size(N,2)
    error('Not enough data collected!')
end

e = [eps*ones(size(xi,1),1)+xi; eps*ones(size(xi,1),1)-xi];
[N, e] = nontrivial_constraints(N, e);
vert = sortrows(con2vert(N,e));

end