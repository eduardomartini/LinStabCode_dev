clear all
% close all
% clc
addpath('aux_matlab');

%% Analisys Parameters
%problem definition
omega       = 0*2*pi;
FDorder = 4;
mm = 25 ;
Re = 30e3 ; 
m  = 20  ;

Ma  = 0.01;
Pr = 0.7;
kappa = 1.4;

% Arnoldi settings
noEigs          = 10;

xRange = [-1.5,-.5] ; nx =  100;
yRange = [0,1]      ; ny =  50;

FDorder=4;
useSimmetry=true;
mesh = SquareMesh(xRange,yRange,nx,ny,FDorder,useSimmetry,false,0);

X=mesh.X; dx=X(2,2)-X(1,1);
Y=mesh.Y; dy=Y(2,2)-Y(1,1);

%% Mesh transform
%parabolic coordinate transformation
Z = (mesh.X+1i*mesh.Y); 
X = -real(Z.^2);
Y = -imag(Z.^2);

disp('Deforming mesh...');
[defmesh] = DeformMesh(mesh,X,Y);

%plot original and physical space meshes
    figure
    subplot(131)
        contourf(mesh.X,mesh.Y,mesh.DW.W); hold on
        plot(mesh.X,mesh.Y,'.k');
        xlabel('$\sigma$');
        ylabel('$\tau$');
        daspect([1,1,1]);
        title('Original mesh and Int. Weights')
           
    subplot(122)
        contourf(defmesh.X,defmesh.Y,defmesh.DW.W); hold on
        plot(defmesh.X,defmesh.Y,'.k');
        xlabel('$x$');
        ylabel('$y$');
        daspect([1,1,1])
        title('Deformerd mesh and Int. Weights')

mesh=defmesh;

%% Set base flow
disp('Setting baseflow...');
BF = ReadBaseFlow(mesh,Re,m,Ma,Pr);

%% Config sponge
%get distance from edges
dist1 = GetMinDistance(mesh,mesh.usedInd([mesh.indexes.left(:)]));
dist2 = GetMinDistance(mesh,mesh.usedInd([mesh.indexes.top(:)]));

%Define sponge strengh, thickness and transition length
spg_str = 15 ; thick1=.15; trans1=1e-2;

%Defines (infinitely smooth) transition function 
f_inf   = @(x) tanh(tan(x)).*(x>-pi/2).*(x<pi/2)+(x>=pi/2)-(x<=-pi/2);
f_inf01 = @(x) (1+f_inf(x*pi/2))/2;

sponge = (1-f_inf01((dist1-thick1)/trans1-1))*spg_str + ...
         (1-f_inf01((dist2-thick1)/trans1-1))*spg_str ;

%% Get linenar operator and set Boundary Conditions
LinProb = GetLinProblem(mesh,BF,sponge,'axy');

DirBCind=unique([mesh.indexes.top;
                 mesh.indexes.left;
                 mesh.indexes.right]);
             
figure;
    plot(mesh.X(:),mesh.Y(:),'.k', ...
    mesh.X(mesh.usedInd(DirBCind)),mesh.Y(mesh.usedInd(DirBCind)),'or');   
    title('Boundary nodes');
    daspect([1,1,1]);
    xlabel('$x$');
    xlabel('$y$');
    
%apply Dir. B.C. for u,v,w, and T 
DirBCind=DirBCind+[1,2,3,4]*mesh.ngp;
DirBCind=DirBCind(:);

LinProb        = SetDirBC(LinProb,DirBCind);

%% Energy Norm
[Wen,invWen] = GetXuEnergyNorm(mesh,BF,'axy');

%% Solve problem

nDOFs   = mesh.ngp*5;
B       = speye(nDOFs);
C       = B;

[G,response,force]=ResAnalysis(LinProb,Wen,invWen,omega,noEigs,B,C);

%% Plot modes and gains

%plot Gains    
f=figure; f.Position (3:4)=[600,600];
    plot(G)
    xlabel('#mode');
    ylabel('$\sigma$');

%plot forcing and response modes
for imode=1:5
    f=figure; f.Position(3:4)=[400,800];
    f.Name=sprintf('Re=%.0fk, Mode #%.0f',Re,imode);
    comps={ LinProb.indexes.u_j,'U';
            LinProb.indexes.v_j,'V';
            LinProb.indexes.w_j,'W';
            LinProb.indexes.rho_j,'RHO';
            LinProb.indexes.T_j,'T';
            };
    for icomp=1:5

    subplot(5,2,1+(icomp-1)*2)
        Uplot = X*nan;
        Uplot(mesh.usedInd)=response(comps{icomp,1},imode);            
        
        contourf(mesh.X,mesh.Y,real(Uplot),64,'linecolor','none');
        title(sprintf('Response %s, w=%.2f', comps{icomp,2},omega))
        colorbar;
        daspect([1,1,1])

    subplot(5,2,2+(icomp-1)*2)
        Vplot = mesh.X*nan;
        Vplot(mesh.usedInd)=force(comps{icomp,1},imode);
        
        contourf(mesh.X,mesh.Y,real(Vplot),64,'linecolor','none');
        title(sprintf('Force %s, w=%.2f', comps{icomp,2},omega))
        colorbar;
        daspect([1,1,1])
    end
end


