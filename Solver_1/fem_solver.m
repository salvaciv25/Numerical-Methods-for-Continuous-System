%%% we want to solve the Convection Diffusion equation
clear all; close all;clc;

%% Define our mesh 
% calling the file CD_mesh
meshdir = 'mesh';

[N_elem, Nodes, NDir, triang,  DirNod, DirVal, xy ]=input_var(meshdir);

TRI = delaunay(xy(:,1),xy(:,2));
trimesh(TRI,xy(:,1),xy(:,2));

%% Computing the local basis

[Area,a,b,c]=LocalB(triang,xy,N_elem);

%% Building the matrix SYSMAT:matrix H, transport matrix B, SUD matrix S
% H: stiff_mat
% B: transpmat
% S: stab_mat
beta = [0 0];
Diff = [8 4];
tau = 0;
[stiff_mat, transpmat, stab_mat]=stiff(N_elem,Nodes,triang,b,c,Area,beta,Diff,tau);

%% solving the linear system by using the Dirichlet boundary condition
SYSMAT= stiff_mat + transpmat + stab_mat;
RHS = zeros(Nodes,1);

%% Imposing the boundary condition
[SYSMAT,RHS]=lifting(SYSMAT,RHS,DirNod,DirVal);
% Parametri per il GMRES
restart = 10; tol = 1e-9; maxit = 20;
setup.type = 'nofill';  % Precondizionamento ILU
setup.milu = 'off';
% Precondizionamento ILU
[L, U] = ilu(SYSMAT, setup);

% Risoluzione del sistema con GMRES
u_h= gmres(SYSMAT, RHS, restart, tol, maxit, L, U);

figure(1)
trisurf(TRI,xy(:,1),xy(:,2),u_h);
xlabel('x');
ylabel('y');
zlabel('u');
hold on
az = -220; 
el = 9;  
view(az, el);
hold off



%% Without stabilization and constant velocity beta=[1 3]
beta = [1 3];
tau = 0;

Diff = [1e-1; 5e-2; 1e-2; 9e-3];


height = floor(sqrt(size(Diff, 1)));
base = ceil(size(Diff, 1) / height);

figure(2)

for i = 1:size(Diff, 1)

    
    [stiff_mat, transpmat, stab_mat,h(i)]=stiff(N_elem,Nodes,triang,b,c,Area,beta,[Diff(i) Diff(i)],tau);

    %% solving the linear system by using the Dirichlet boundary condition
    SYSMAT= stiff_mat + transpmat + stab_mat;
    RHS = zeros(Nodes,1);

    %% Imposing the boundary condition
    [SYSMAT,RHS]=lifting(SYSMAT,RHS,DirNod,DirVal);
    % Parametri per il GMRES
    restart = 10; tol = 1e-9; maxit = 20;
    setup.type = 'nofill';  % Precondizionamento ILU
    setup.milu = 'off';
    % Precondizionamento ILU
    [L, U] = ilu(SYSMAT, setup);

    % Risoluzione del sistema con GMRES
    u_h= gmres(SYSMAT, RHS, restart, tol, maxit, L, U);


    subplot(height, base, i)
    trisurf(TRI,xy(:,1),xy(:,2),u_h, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('u');
    title(sprintf('K=%g', Diff(i)));

    az = -220; 
    el = 9; 
    view(az, el);
end


%% with stabilization 
Diff = 1e-2;
beta = [1, 3];
tau = [1e-2; 0.1; 0.5; 1];

height = floor(sqrt(size(tau, 1)));
base = ceil(size(tau, 1) / height);

figure(3)

for i = 1:size(tau, 1)

    
    [stiff_mat, transpmat, stab_mat]=stiff(N_elem,Nodes,triang,b,c,Area,beta,[Diff Diff],tau(i));

    %% solving the linear system by using the Dirichlet boundary condition
    SYSMAT= stiff_mat + transpmat + stab_mat;
    RHS = zeros(Nodes,1);

    %% Imposing the boundary condition
    [SYSMAT,RHS]=lifting(SYSMAT,RHS,DirNod,DirVal);
    % Parametri per il GMRES
    restart = 10; tol = 1e-9; maxit = 20;
    setup.type = 'nofill';  % Precondizionamento ILU
    setup.milu = 'off';
    % Precondizionamento ILU
    [L, U] = ilu(SYSMAT, setup);

    % Risoluzione del sistema con GMRES
    u_h= gmres(SYSMAT, RHS, restart, tol, maxit, L, U);
    
    subplot(height, base, i)
    trisurf(TRI,xy(:,1),xy(:,2),u_h, 'FaceColor', 'interp', 'EdgeColor', 'black');
    xlabel('x');
    ylabel('y');
    zlabel('u');
    title(sprintf('tau=%i'), tau(i));

    az = -220; 
    el = 9;  
    view(az, el);
end