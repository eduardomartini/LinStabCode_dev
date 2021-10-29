clear all
% close all
% clc
% restoredefaultpath
addpath('aux_matlab');

%% Analisys Parameters
%problem definition
omega       = 0*2*pi;
FDorder = 4;
Re = 10e3 ; 
alpha  = 1  ;

Ma  = 0.01;
Pr = 0.7;
kappa = 1.4;

% Arnoldi settings
noEigs          = 10;

xLengh = .01 ; nx = 1;
yLenght = 1  ; ny =200;

% FD
FDorder=4;
useSimmetry=false;
periodic=true;
xRange=[0,xLengh];
yRange=[-yLenght,yLenght];

mesh = CreateRectangularPeriodicMesh(xRange,yRange,nx,ny,FDorder,useSimmetry);

X=mesh.X; 
if nx>1
    dx=X(1,2)-X(1,1);
else
    dx=1;
end
Y=mesh.Y; dy=abs(Y(2,1)-Y(1,1));
% 
indexes.bottom    = find(Y(mesh.usedInd)<min(Y(:))+dy/2);
indexes.top       = find(Y(mesh.usedInd)>max(Y(:))-dy/2);
indexes.left      = find(X(mesh.usedInd)<min(X(:))+dx/2);
indexes.right     = find(X(mesh.usedInd)>max(X(:))-dx/2);
 
%% Set base flow
disp('Setting baseflow...');

BF.Ma= Ma;
BF.Pr= Pr;

BF.RHO      = X*0+1;
% BF.U        = 0.5*(tanh(Y)+1);
BF.U        = 1-Y.^2;
BF.V        = X*0;
BF.W        = X*0;

BF.T        = X*0+1;
BF.MU       = X*0+1;
BF.kappa    = kappa;
BF.cv       = 1/(BF.kappa*(BF.kappa-1)*Ma^2);
BF.c1       = (BF.kappa-1)*Re*Pr*Ma^2;
BF.c2       = BF.kappa*Ma^2;


BF.dMUdT = X*0;
BF.dmudT = X*0;

BF.d2MUdT2 = X*0;
BF.d2mudT2 = X*0;

sponge=zeros(size(mesh.X));

%% Get linenar operator and set Boundary Conditions
[LHS,RHS,compIndex] = GetLHSRHS_axy(mesh,BF,alpha,Re,sponge,'2D');

DirBCind=unique([indexes.bottom;
                 indexes.top]);
DirBCind=DirBCind+[1,2,3,4]*mesh.ngp;
DirBCind=DirBCind(:);

[LHS,RHS] = SetDirichletBC(LHS,RHS,DirBCind);
L0    = -(RHS/1i)\LHS;

mem = whos('L0');
fprintf('Memory for L0 storage : %0.4f Gb\n', mem.bytes/1e9);

%% Eigenmodes
tic()
lambda_tar = (0)*alpha;

[V,lambda]=eigs(L0/-1i,200,lambda_tar);

fprintf("Eigs took : %.1f s\n",toc());
lambda=diag(lambda);

[~,order]=sort(abs(lambda-lambda_tar),'ascend');
lambda=lambda(order);
V=V(:,order);

figure
    plot(lambda,'o');hold on;
    xlim([-1,3]);
    ylim([-1,.3]);
    legend