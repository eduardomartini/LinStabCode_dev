%set memory log
% profile clear
% profile -memory on

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex')
clear variables
% close all
clc
addpath('aux_matlab');

% Physical parameters
baseFlow.Re     = 5e2;          % Reynolds number
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
verbose     = false;         % visualize grid, base flow and results

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
spongeAmp   = 15;
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


%% Eigenvalues 
neigs = 1;


opts.tol    = eps;
opts.disp   = 2;
opts.issym  = false;
opts.isreal = false;

% using standard approaches
L1 = L0;
L1(idx_dirchlet,:)=[];
L1(:,idx_dirchlet)=[];
lambdatar = 15;
tic
    [V,lambda_std] = eigs(L1,neigs,lambdatar, opts);
%     lambda_std = zeros(neigs)
time_eig_std  = toc;
disp(['Elapsed time for classical approach, ' num2str(time_eig_std) 's'])
lambda_std=diag(lambda_std);



% using time marching

    % mask matrix to impose constrains
H = speye(size(L0));
H(idx_dirchlet,idx_dirchlet)=0;


    % setup scheme
dt = 5e-1;

maxIter = 1;
n_block        = [maxIter];
i_block_output = [0];

tic
%     H = speye(size(L1));

%     [TM_setup,TM_setup_adj] = TM_EulerImplicit_Setup(H,L1,dt,[],[],1e-4,0);
    [TM_setup,TM_setup_adj] = TM_EulerImplicit_Setup(H,L0,dt,[],[],1e-4,0);
%     TM_setup = TM_BDF2_Setup(H,L0,dt,[],[],1e-4,2);
    time_eig_TM_part1 = toc();
    Dir = @(q0) TM(q0,@(t) zeros(size(H,1),1) ,maxIter,TM_setup   ,i_block_output,n_block,1e-6);
    Adj = @(q0) TM(q0,@(t) zeros(size(H,1),1) ,maxIter,TM_setup_adj,i_block_output,n_block,1e-6);

%     ritz1 = lambda_std;
   [V1,ritz1] = eigs(Dir, size(H,1), neigs, 'lm', opts);
   ritz1=diag(ritz1);
   ritz2 = eigs(Adj, size(H,1), neigs, 'lm', opts);
time_eig_TM = toc;
lambda_tm1 = (ritz1-1)./(dt*ritz1);
lambda_tm2 = (ritz2-1)./(dt*ritz2);
disp(['Elapsed time for time marching approach, ' num2str(time_eig_TM) 's'])

disp([size(L0,1)])
disp([lambda_std,lambda_tm1,lambda_tm2])

%%
    ndofs = size(L0,1);
    used_dofs = 1:ndofs;  used_dofs(idx_dirchlet)=[];
    U = zeros(ndofs,1);
    U=V1;
    figure('name','Resolvent forcing and response modes')
    vars = {real(U(idx.rho_j)) ,'$\rho$'; 
            real(U(idx.u_j  )) ,'$u$'; 
            real(U(idx.v_j  )) ,'$v$'; 
            real(U(idx.w_j  )) ,'$w$'; 
            real(U(idx.T_j  )) ,'$T$' };
    plotFlow(mesh.X,mesh.Y,vars,3,2)


fprintf('Time (std/TM) : %3.3e / %3.3e / %3.3e , error : %3.1e\n',time_eig_std,time_eig_TM_part1,time_eig_TM, abs(sort(lambda_std)-sort(lambda_tm)))