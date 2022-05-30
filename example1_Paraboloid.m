set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex')
clear variables
close all
clc
addpath('aux_matlab');

% Physical parameters
baseFlow.Re     = 5e3;          % Reynolds number
baseFlow.Ma     = 0.1;          % Mach number
baseFlow.Pr     = 0.7;          % Prandtl number
baseFlow.kappa 	= 1.4;          % heat capacity ratio
baseFlow.T_0 	= 293.15;       % temperature

% Perturbation & EVP parameters
freq        = 0;            % frequency
omega       = freq*2*pi;    % angular frequency
m           = 30;           % azimuthal wave number
nEig        = 3;            % Arnoldi method number of eigenvalues

% Domain & grid
Nr          = 50;           % # of grid points (radial)
Nz          = 50;           % # of grid points (streamwise)
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
xrange          = [-1 0 ];  % domain range in x
yrange          = [ 0 1 ];  % domain range in y

cmesh           = CreateMesh(xrange,yrange,Nz,Nr,FDorder, ...     
                             y_symmetry,x_periodicity,alpha); %construct mesh
                     
x   = cmesh.X;           % x,y: Cartesian grid coordinates
y   = cmesh.Y;

% Grid transformation to parabolic sement in physical domain
d   = 0.2;              % wall-normal thickness parameter
x   = x*d-0.5;
z   = (x+1i*y);
X   = -real(z.^2);      % X,Y: Parabolic grid coordinates
Y   = -imag(z.^2);

%
mesh    = DeformMesh(cmesh,X,Y);
mesh.W  = mesh.W.*mesh.Y;   % cylindrical volume element (dV=int(rdrdz))

% Plot meshes in comutational and physical domains
if verbose
    figure
    
    subplot(1,2,1)
    title('Computatinal domain')
    hold on
    plot(cmesh.X, cmesh.Y, '-k', cmesh.X', cmesh.Y', '-k','HandleVisibility','off');
    for bondaries= fields(mesh.idx)'
        ids = mesh.idx.(bondaries{1});
        plot(cmesh.X(ids), cmesh.Y(ids),'linewidth',2);
    end
    xlabel('$\xi$');
    ylabel('$\eta$');
    legend(fields(mesh.idx))
    axis equal tight
    
    subplot(1,2,2)
    title('Physical domain')
    hold on
    plot(mesh.X, mesh.Y, '-k', mesh.X', mesh.Y', '-k','HandleVisibility','off');
    for bondaries= fields(mesh.idx)'
        ids = mesh.idx.(bondaries{1});
        plot(mesh.X(ids), mesh.Y(ids),'linewidth',2);
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
baseFlow    = example_1_readbaseflow(mesh,baseFlow); % Custumize the function to read your baseflow

% Viscosity (via Sutherland) and heat conductivity (via constant Prandtl number)
[baseFlow]  = sutherland_air(baseFlow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sponge                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sponge defined in the origial, non-deformed, mesh.
xs          = -0.5-d*4/5;
xs_trans    = d/10;
ys          = 0.9;
ys_trans    = 0.01;
spongeAmp   = 2;
mesh.sponge = spongeAmp/2*max(tanh( (y-ys)/ys_trans)+1,  ...
    tanh(-(x-xs)/xs_trans)+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize base flow                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    figure('name','Base Flow')
    vars = {baseFlow.RHO ,'$\rho$';     baseFlow.U,'$U$';
            baseFlow.V,'$V$';           baseFlow.W,'$W$';
            baseFlow.T,'$T$';           mesh.sponge,'sponge';
            baseFlow.MU,'$\mu$';        baseFlow.dmudT,'$\frac{d\mu}{dT}$';
            baseFlow.d2mudT2                                ,'$\frac{d^2\mu}{dT^2}$';
            reshape(mesh.Dy*(baseFlow.W(:)),Nr,Nz)          ,'$\frac{dW}{dy}$';
            reshape(mesh.Dx*(baseFlow.W(:)),Nr,Nz)         ,'$\frac{dW}{dx}$';
            reshape(mesh.Dx*(mesh.Dy*baseFlow.W(:)),Nr,Nz) ,'$\frac{d^2W}{dydx}$'};
    plotFlow(mesh.X,mesh.Y,vars,4,3)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up linear operator incl. boundary conditions                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dqdt = L0 q
%%
[L0,idx] 	= GetLinProblem(mesh,baseFlow,'axy',m);

% Enforce Dirichlet b,c on the top, right and left boundaries, for u,v,w
% and T
borders='ltr';  vars = 'uvwT';
[L0,idx_dirchlet] = BC_Dirichlet(L0,idx,borders,vars);

[W,invW] = GetXuEnergyNorm(mesh,baseFlow,'axy');

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
    vars = {real(V(idx.u_j,1)) ,'$f_u^{(1)}$'; real(U(idx.u_j,1)) ,'$u^{(1)}$';real(V(idx.u_j,2)) ,'$f_u^{(2)}$'; real(U(idx.u_j,2)) ,'$u^{(2)}$';real(V(idx.u_j,3)) ,'$f_u^{(3)}$'; real(U(idx.u_j,3)) ,'$u^{(3)}$';};
    plotFlow(mesh.X,mesh.Y,vars,3,2)
end
