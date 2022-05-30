set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex')
clear variables
close all
clc
addpath('aux_matlab');

% Physical parameters
baseFlow.Re     = 100;          % Reynolds number
baseFlow.Ma     = 0.01;         % Mach number
baseFlow.Pr     = 0.7;          % Prandtl number
baseFlow.kappa 	= 1.4;          % heat capacity ratio
baseFlow.T_0 	= 293.15;       % temperature

% Perturbation & EVP parameters
freq        = 0;            % frequency
omega       = freq*2*pi;    % angular frequency
beta        = 0;            % azimuthal wave number
nEig        = 3;            % Arnoldi method number of eigenvalues

% Domain & grid
Nx          = 30*2+1;           % # of grid points (radial)
Ny          = 40*3+4;           % # of grid points (streamwise)
FDorder     = 4;            % finite difference order of accuracy


% Flags
verbose     = true;         % visualize grid, base flow and results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create mesh and obtain differentiation matrices                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cartesian mesh in computational domain
y_symmetry      = true;     % use symmetry on y coordinate around y=0 (for axysymmetric problems)
x_periodicity   = false;    % use periodic b.c. on x
alpha           = .0;       % spatial filter coefficient
xrange          = [ -1 2 ];  % domain range in x
yrange          = [ -1 1 ];  % domain range in y

cmesh           = CreateMesh(xrange,yrange,Ny,Nx,FDorder, ...     
                             y_symmetry,x_periodicity,alpha); %construct mesh
                     
x   = cmesh.X;           % x,y: Cartesian grid coordinates
y   = cmesh.Y;


mask = (y<0) & ( (x<0) | (x>1) ) ;  

mesh = MeshMask(cmesh,mask);

%%
x = mesh.X;
y = mesh.Y;
p = mesh.usedInd;

%%

% Plot meshes in comutational and physical domains
if verbose
    figure
   
    title('Physical domain')
    hold on
    ind = mesh.usedInd;
    plot(mesh.X(ind), mesh.Y(ind), '.k','HandleVisibility','off');
    for bondaries= fields(mesh.idx)'
        ids = mesh.idx.(bondaries{1});
        plot(mesh.X(p(ids)), mesh.Y(p(ids)),'o');
    end
    legend(fields(mesh.idx),'Location', 'Best')
    xlabel('$x$');
    ylabel('$y$');
    axis equal tight
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Base flow import                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import base flow

disp('Setting baseflow...');
baseFlow    = ReadStepBaseFlow(mesh,baseFlow); % Custumize the function to read your baseflow

% Viscosity (via Sutherland) and heat conductivity (via constant Prandtl number)
%[baseFlow]  = sutherland_air(baseFlow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sponge                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sponge defined in the origial, non-deformed, mesh.
xs_trans    = .1;
ys_trans    = .1;
spongeAmp   = 2;
ind = mesh.usedInd;

mesh.sponge = nan(size(mesh.X));
mesh.sponge(ind) = ...
        spongeAmp/2*max(tanh( ( y(ind)-max(y(:)) )/ys_trans)+1, ...
                        tanh(-( x(ind)-min(x(:)) )/xs_trans)+1 + ...
                        tanh( ( x(ind)-max(x(:)) )/xs_trans)+1  );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize base flow                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    figure('name','Base Flow')
    vars = {baseFlow.RHO ,'$\rho$';     baseFlow.U,'$U$';
            baseFlow.V,'$V$';           baseFlow.W,'$W$';
            baseFlow.T,'$T$';           mesh.sponge,'sponge';
            baseFlow.MU,'$\mu$'; 
            mesh.W,'$W$'; 
            
         %   reshape(mesh.Dy*(baseFlow.W(:)),Nx,Ny)          ,'$\frac{dW}{dy}$';
         %   reshape(mesh.Dx*(baseFlow.W(:)),Nx,Ny)         ,'$\frac{dW}{dx}$';
         %   reshape(mesh.Dx*(mesh.Dy*baseFlow.W(:)),Nx,Ny) ,'$\frac{d^2W}{dydx}$'
         };
    plotFlow(mesh.X,mesh.Y,vars,4,3)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up linear operator incl. boundary conditions                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dqdt = L0 q
%%
[L0,idx] 	= GetLinProblem(mesh,baseFlow,'2D',beta);

% Enforce Dirichlet b,c on the top, right and left boundaries, for u,v,w
% and T
borders='lbtrm';  vars = 'uvwT';
[L0,idx_dirchlet] = BC_Dirichlet(L0,idx,borders,vars);

[W,invW] = GetXuEnergyNorm(mesh,baseFlow,'2D');
W = 1;
invW=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Resolvent analysis                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dq/dt = L*q + B*v
% u = C*q
% v: input/forcing, u:output/response
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input/forcing matrix
B                       = ones(mesh.ngp*5,1);
B([idx.T_j;idx.rho_j])  = 0;	% no temperature and density forcing (momentum only)
B( idx_dirchlet  )      = 0;    % no forcing on dirichlet b.c. 
B                       = spdiags(B,0,mesh.ngp*5,mesh.ngp*5);

% Output/observable matrix
C                       = B;    % ignore output temperature, density and dir.b.c. (same as for inputs)

% Compute optimal forcings and responses
[S,U,V]                 = resolvent(L0,W,invW,omega,nEig,B,C,mesh.filters);

%% Plot modes and gains
if verbose
    figure('name','Mode gains')
    bar(S.^2);
    xlabel('mode');
    ylabel('gain');
    title('Resolvent gains')
    
    figure('name','Resolvent forcing and response modes')
    vars = {real(V(idx.u_j,1)) ,'$f_u^{(1)}$'; 
            real(U(idx.u_j,1)) ,'$u^{(1)}$';
            real(V(idx.u_j,2)) ,'$f_u^{(2)}$'; 
            real(U(idx.u_j,2)) ,'$u^{(2)}$';
            real(V(idx.u_j,3)) ,'$f_u^{(3)}$'; 
            real(U(idx.u_j,3)) ,'$u^{(3)}$';};
    plotFlow(mesh.X,mesh.Y,vars,3,2,mesh.usedInd)
end
