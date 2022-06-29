set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex'); set(groot,'defaultTextInterpreter','latex')
clear 
% close all
clc
addpath('aux_matlab');

measure_memory = false ;


if measure_memory
    set memory log
    profile clear
    profile -memory on
end

% Physical parameters
baseFlow.Re     = 28000;          % Reynolds number
baseFlow.Ma     = 6;          % Mach number
baseFlow.Pr     = 0.7;          % Prandtl number
baseFlow.kappa 	= 1.4;          % heat capacity ratio
baseFlow.T_0 	= 273.15 ;       % temperature

% Perturbation & EVP parameters
freq        = 0;            % frequency
omega       = freq*2*pi;    % angular frequency for spatial stability
alpha       = .6283;        % stream-wize wavenumber for temporal stability
nEig        = 15;            % Arnoldi method number of eigenvalues

% Domain & grid
Nx          = 50;           % # of grid points (radial)
Ny          = 25;           % # of grid points (streamwise)
FDorder     = 4;            % finite difference order of accuracy


% Flags
verbose     = true;         % visualize grid, base flow and results

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create mesh and obtain differentiation matrices                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cartesian mesh in computational domain
y_symmetry      = false;         % use symmetry on y coordinate around y=0 (for axysymmetric problems)
x_periodicity   = true;          % use periodic b.c. on x
alpha_filter    = .0;            % spatial filter coefficient
xrange          = [0  pi ];  % domain range in x
yrange          = [ -1   1 ];  % domain range in y

cmesh           = CreateMesh(xrange,yrange,Nx,Ny,FDorder, ...     
                             y_symmetry,x_periodicity,alpha_filter); %construct mesh
                     
X   = cmesh.X;           % x,y: Cartesian grid coordinates
Y   = cmesh.Y;

yi= 2.0; %dns_scale*delta ; % Hanifi ?
ymax=8 ; % Hanifi used from 100 to 500
a= yi*ymax/(ymax-2*yi);
b= 1+2*a/ymax;

Y=a*(1+Y)./(b-Y);  % Stretch to semi-infinite domain

%
mesh    = DeformMesh(cmesh,[],Y,true);

% mesh    = DeformMesh(cmesh,X,Y,true);
% mesh    = cmesh;

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
baseFlow    = example_5_readbaseflow(mesh,baseFlow); % Custumize the function to read your baseflow

% Viscosity (via Sutherland) and heat conductivity (via constant Prandtl number)
% [baseFlow]  = sutherland_air(baseFlow);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sponge                                                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sponge defined in the origial, non-deformed, mesh.
ys_trans    = 1;
ys          = ymax-2;
spongeAmp   = .1;
mesh.sponge = spongeAmp/2*(tanh( (Y-ys)/ys_trans)+1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize base flow                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    figure('name','Base Flow')
    vars = {baseFlow.RHO ,'$\rho$';     baseFlow.U,'$U$';
            baseFlow.V,'$V$';           baseFlow.W,'$W$';
            baseFlow.T,'$T$';           mesh.sponge,'sponge';
            baseFlow.MU,'$\mu$';        baseFlow.dMUdT,'$\frac{d\mu}{dT}$';
            baseFlow.d2MUdT2                                ,'$\frac{d^2\mu}{dT^2}$';
            reshape(mesh.Dy*(baseFlow.W(:)),Nx,Ny)          ,'$\frac{dW}{dy}$';
            reshape(mesh.Dx*(baseFlow.W(:)),Nx,Ny)         ,'$\frac{dW}{dx}$';
            reshape(mesh.Dx*(mesh.Dy*baseFlow.W(:)),Nx,Ny) ,'$\frac{d^2W}{dydx}$'};
    plotFlow(mesh.X,mesh.Y,vars,4,3)
    drawnow
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set up linear operator incl. boundary conditions                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    dqdt = L0 q
%%
[L0,idx] 	= GetLinProblem(mesh,baseFlow,'2D',alpha);
% Enforce Dirichlet b,c on the top, right and left boundaries, for u,v,w
% and T
borders='tb';  vars = 'uvwT';
[L0,idx_dirchlet] = BC_Dirichlet(L0,idx,borders,vars);


%%     Compute eigenvalues       %%
omega_tar   = 0.5498+1.1i;   % value around which the eigenvalues will be seek

%-------- Classical approach ------

% Remove constrains lines/rows and get L1 : \omega \hat q = L1 \hat q
L1 = L0/-1i;
L1(idx_dirchlet,:)=0;
L1(:,idx_dirchlet)=0;
L1(idx_dirchlet,idx_dirchlet) = -1e4;

%setup eigs options
opts.tol    = 1e-8;
opts.disp   = 2;
opts.issym  = false;
opts.isreal = false;
opts.p      = 150;


% compute eigs
tic;
[V,omega] = eigs(L1,nEig,omega_tar,opts);
omega=diag(omega);
time_eigsMatlab = toc;

% compute eigs using predifined function
for type={'ilu'}
    disp(type{1})
    optsInv.type=type{1};
    optsInv.verbose=false;
    optsInv.tol       =1e-4  ; 
    optsInv.toliLU    =1e-6  ; 
 
    tic
    [invA_fun,invA_T_fun] = GetInverseFunction(L1-speye(size(L1))*omega_tar,optsInv);
%     lambda2.(type{1})=nan;
%     [V2,lambda2.(type{1})] = eigs(invA_fun,size(L1,1),nEig,omega_tar,opts); 
    lambda2.(type{1})=diag(lambda2.(type{1}));
    time.(type{1}) = toc;
end

% compute eigs using a block-arnoldi method with parallelization 
%       ( still being optimized )

tic
n = 100; %opts.p;
m = n/2;
p = 8;
% [invA_fun,invA_T_fun] = GetInverseFunction(L1-speye(size(L1))*omega_tar);
[V3,ritz3] = block_eigs2(invA_fun, size(L1,1), nEig,n,m,p,opts.tol);
% ritz3=nan;
time_eigsblock = toc;
lambda3 = 1./(ritz3)-omega_tar;   %revert shift and invert 
%%
    figure('name','Temporal Spectra');
        plot(real(omega),imag(omega)*0,'ro',real(lambda2.ilu),imag(lambda2.ilu),'b+',real(lambda3),imag(lambda3)*0,'g*')


%% Plot and compare eigenvalues. Show some eigen vectors
if verbose    
    [~,iplot] = sort(imag(omega),"descend");
    iplot=iplot(1:4);
    
%     [~,iplot]=min(abs(lambda-(.9*alpha+.10i)));
    figure('name','Temporal Spectra');
        plot(real(omega),imag(omega),'ro',real(lambda2),imag(lambda2),'b+',real(lambda3),imag(lambda3),'g*')
        legend( ['eigs (' num2str(time_eigsMatlab,'%3.1f') 's)'] , ... 
                ['eigs (pre.def. fun.) (' num2str(time_eigsMatlab_fun,'%3.1f') 's)'] , ...
                ['block-arnoldi (' num2str(time_eigsblock,'%3.1f') 's)'])

%         plot(real(lambda)/alpha,imag(lambda),'ob');
        hold on
        for iiplot=iplot'
            plot(real(omega(iiplot)),imag(omega(iiplot)),'o','handlevisibility','off');
        end
        grid on
        xlabel('$\omega_r$');
        ylabel('$\omega_i$');
    
    ndofs = size(L0,1);
    used_dofs = 1:ndofs;  used_dofs(idx_dirchlet)=[];
    U = zeros(ndofs,1);
    for iiplot=iplot'
        U=V(:,iiplot);
        figure('name','Eigenmode')
        vars = {real(U(idx.rho_j)) ,'$\rho$'; 
                real(U(idx.u_j  )) ,'$u$'; 
                real(U(idx.v_j  )) ,'$v$'; 
                real(U(idx.w_j  )) ,'$w$'; 
                real(U(idx.T_j  )) ,'$T$' };
        
        axs = plotFlow(mesh.X,mesh.Y,vars,3,2,[],101,'linecolor','none');
        for i=1:length(axs)
            axs(i).YLim=[0,4];
        end
    end
end
    
    
%% Spatial Stability Analysis
[L,Lw,R,idx,l0,r0,r1,r2] 	= GetSpatialLinProblem(mesh,baseFlow,'2D');

% Enforce Dirichlet b,c on the top, right and left boundaries, for u,v,w
% and T
borders='tb';  vars = 'uvwT';
[~,idx_dirchlet] = BC_Dirichlet(L,idx,borders,vars);

w = 2*pi*0.09; 
alpha_target = alpha;

LL = L+Lw*w ;
RR = R      ;
LL(idx_dirchlet,:) = [];
LL(:,idx_dirchlet) = [];

RR(idx_dirchlet,:) = [];
RR(:,idx_dirchlet) = [];

alphas = eigs(LL,RR,nEig,alpha_target,opts);

%%
figure('name','Spatial spectra')
    plot(real(alphas),imag(alphas),'ob')
    xlabel('$\alpha_r$');
    ylabel('$\alpha_i$');
    grid on;
    
abort
%% Compute optimal forcings and responses %%
%% Classical approach
tic
[S,U,V]                 = resolvent(L0,W,invW,omega,nEig,B,C,mesh.filters);
time_std = toc;
%% Time marching methods 
% 
H = speye(size(L0));
H(idx_dirchlet,idx_dirchlet)=0;
dt=5e-1;

tic
[TM_setup,TM_setup_adj] = TM_EulerImplicit_Setup(H,L0,dt,[],[],1e-4,2);
time_iLU = toc;

tic
nIter = 5;
tol=1e-3;
deltaF = 1;
nfreqs =1 ;
mRSVD=3;
[Stm,~,fList,SS_conv] = TM_Resolvent(TM_setup,TM_setup_adj,deltaF,nfreqs,nIter,tol,W,invW,B,C,mRSVD);
time_tm = toc;

figure('name','Gains convergence')
    plot(1:nIter, SS_conv.',1:nIter,repmat(S,1,nIter),'k:') 
    xlabel('iterations');
    ylabel('gain');
    title('Gain convergence')
    
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

%% Print memory log
if measure_memory
    profile report
    p = profile('info');

    %%
    n= length(p.FunctionTable);
    for i=1:n
        if strcmp(p.FunctionTable(i).FunctionName , 'resolvent')
            MatrixTime = p.FunctionTable(i).TotalTime;
            MatrixMem  = p.FunctionTable(i).PeakMem;        
        end
        if strcmp(p.FunctionTable(i).FunctionName , 'TM_Resolvent')
            MatrixfreeTime = p.FunctionTable(i).TotalTime;
            MatrixfreeMem  = p.FunctionTable(i).PeakMem;        
        end
    end
    dlmwrite('example5_mem_history.csv',[Nr,Nz,MatrixTime,MatrixMem,MatrixfreeTime,MatrixfreeMem,dt],'-append')
end
    