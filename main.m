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

xRange = [-1.5,-.5] ; nx =  120;
yRange = [0,1]     ; ny =  30;

FDorder=4;
useSimmetry=true;
mesh = SquareMesh(xRange,yRange,nx,ny,FDorder,useSimmetry);

X=mesh.X; dx=X(2,2)-X(1,1);
Y=mesh.Y; dy=Y(2,2)-Y(1,1);
% 

mask =  (mesh.X<-1 & mesh.X<-.6 & mesh.Y>.75 ) | ...
       (mesh.X>-1.3 & mesh.X<-1.1 & mesh.Y<.35 )  ;
n=2;
% mask = ( (mesh.X-min(mesh.X(:))) - 1*(mesh.Y-mesh.Y(n+1,1)) < 0  ) ;
mask(1:n,:)=repmat(mask(n+1,:),n,1);
% mask(end+1-(1:n),:)=repmat(mask(end-n,:),n,1);

mask(:)=0; 

[mesh,maskedge] = MeshMask(mesh,mask);

indexes.bottom    = find(Y(mesh.usedInd)<min(Y(:))+dy/2);
indexes.top       = find(Y(mesh.usedInd)>max(Y(:))-dy/2);
indexes.left      = find(X(mesh.usedInd)<min(X(:))+dx/2);
indexes.right     = find(X(mesh.usedInd)>max(X(:))-dx/2);


k=40;
F = sin(k*mesh.Y)/k^2;
dF=-sin(k*mesh.Y);

F(mask)=nan;
dF(mask)=nan;

subplot(211)
    dxF = mesh.X*nan;
    dxF(mesh.usedInd) = mesh.DW.D2y_asymm*F(mesh.usedInd);
    contourf(mesh.X,mesh.Y,dxF-dF);colorbar
subplot(212)
    hold off
    contourf(mesh.X,mesh.Y,F);colorbar;
    hold on;

figure
    plot(mesh.X(mesh.usedInd),mesh.Y(mesh.usedInd),'k.');
    hold on;
    plot(mesh.X(mesh.usedInd(maskedge)),mesh.Y(mesh.usedInd(maskedge)),'ob');
    plot(mesh.X(mesh.usedInd(indexes.top)),mesh.Y(mesh.usedInd(indexes.top)),'or');
    plot(mesh.X(mesh.usedInd(indexes.left)),mesh.Y(mesh.usedInd(indexes.left)),'og');

%%

alpha_filter=0;
[FILTER,FILTER_ct] = GetFilter(mesh,alpha_filter);
         
% FILTER    = @(x) x;
% FILTER_ct = @(x) x;


% Transform mesh
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
        daspect([1,1,1])        
           
    subplot(122)
        contourf(defmesh.X,defmesh.Y,defmesh.DW.W); hold on
        plot(defmesh.X,defmesh.Y,'.k');
        xlabel('$x$');
        ylabel('$y$');
        daspect([1,1,1])

    
mesh=defmesh;
% clear defmesh;

%% Set base flow
disp('Setting baseflow...');

if Re>=1e3
    nekFile=sprintf('bf_meshNL_Re%.0fk.fld',Re/1e3);
else
    nekFile=sprintf('bf_meshNL_Re%.1fk.fld',Re/1e3);
end
% nekFile='paraboloidFlow_Re5k.fld';

data = readnek(nekFile);
Xnek = data(:,:,1);
Ynek = data(:,:,2);
Unek = data(:,:,3);
Vnek = data(:,:,4);
clear data;

BF.RHO      = X*0+1;
BF.U        = X*0;
BF.V        = X*0;
BF.W        = X*0;

BF.W(:) = griddata(Xnek(:),Ynek(:),Unek(:),X(:),Y(:));
BF.U(:) = griddata(Xnek(:),Ynek(:),Vnek(:),X(:),Y(:));

BF.T        = X*0+1;
BF.MU       = X*0+1;
BF.kappa    = 1.4;
BF.cv       = 1/(BF.kappa*(BF.kappa-1)*Ma^2);
BF.c1       = (BF.kappa-1)*Re*Pr*Ma^2;
BF.c2       = BF.kappa*Ma^2;


BF.dMUdT = X*0;
BF.dmudT = X*0;

BF.d2MUdT2 = X*0;
BF.d2mudT2 = X*0;

%Setup sponge
%sponge parameters

dist1 = GetMinDistance(mesh,mesh.usedInd([indexes.left(:);maskedge(:)]));
dist2 = GetMinDistance(mesh,mesh.usedInd([indexes.top(:)]));

spg_str = 15 ;
thick1=.15; trans1=1e-2;

f_inf   = @(x) tanh(tan(x)).*(x>-pi/2).*(x<pi/2)+(x>=pi/2)-(x<=-pi/2);
f_inf01 = @(x) (1+f_inf(x*pi/2))/2;

sponge_ =  (1-f_inf01((dist1-thick1)/trans1-1))*spg_str + ...
          (1-f_inf01((dist2-thick1)/trans1-1))*spg_str ;
sponge=sponge_*nan;
sponge(mesh.usedInd)=sponge_(mesh.usedInd);

figure
    contourf(mesh.X,mesh.Y,sponge);colorbar
    hold on;
    plot(X(:),Y(:),'.w','Markersize',1)

    daspect([1,1,1]);

%% 
% abort
%% Get linenar operator and set Boundary Conditions
[LHS,RHS,compIndex] = GetLHSRHS_axy(mesh,BF,m,Re,sponge,'axy');

DirBCind=unique([indexes.top;
                 indexes.left;
                 indexes.right;
                 maskedge]);
             
figure;
    plot(mesh.X(:),mesh.Y(:),'.k', ...
         mesh.X(mesh.usedInd(DirBCind)),mesh.Y(mesh.usedInd(DirBCind)),'or');   
     hold on;
     contourf(mesh.X,mesh.Y,sponge)

DirBCind=DirBCind+[1,2,3,4]*mesh.ngp;

DirBCind=DirBCind(:);

[LHS,RHS] = SetDirichletBC(LHS,RHS,DirBCind);


%% Energy Norm
intWeight = mesh.DW.W .*mesh.Y; %int () r dr dz
p=mesh.usedInd;
F  = ...
    0.5*[
    (BF.T(p) ./(BF.RHO(p)*BF.kappa*Ma^2)).* intWeight(p)
     BF.RHO(p).*intWeight(p)
     BF.RHO(p).*intWeight(p)
     BF.RHO(p).*intWeight(p)
    (BF.RHO(p)./(BF.kappa*(BF.kappa-1)*BF.T(p)*Ma.^2)).*intWeight(p)
    ];
invF = 1./F;

%% Solve problem

nDOFs   = mesh.ngp*5;
B       = ones(nDOFs,1);

B(DirBCind)=0;

B=spdiags(B,0,nDOFs,nDOFs);
C=B;

L0    = -(RHS/1i)\LHS + spdiags(repmat(sponge(mesh.usedInd),5,1),0,nDOFs,nDOFs);

gain  = cell(length(omega),1);
U_out = cell(length(omega),1);
V_in  = cell(length(omega),1);

[G,response,force]=ResAnalysis(L0,F,invF,omega,noEigs,B,C,FILTER,FILTER_ct);
%         [G,response,force]=ResAnalysisTimeMarching(L0,1e-1,1e-6,F,invF,omega,noEigs,B,C,DirBCind,1e-6,1e-6);



w=omega;



%% Plot modes and gains

    
f=figure; f.Position (3:4)=[600,600];
    plot(G)
%     title(sprintf('Nr=%.0f, Nz=%0.f',Ny,Nx));
    xlabel('#mode');
    ylabel('$\sigma$');


for imode=1:5
    f=figure; f.Position(3:4)=[400,800];
    f.Name=sprintf('Re=%.0fk, Mode #%.0f',Re,imode);
    comps={ compIndex.u_j,'U';
            compIndex.v_j,'V';
            compIndex.w_j,'W';
            compIndex.rho_j,'RHO';
            compIndex.T_j,'T';
            };
    for icomp=1:5

    subplot(5,2,1+(icomp-1)*2)
        Uplot = X*nan;
        Uplot(mesh.usedInd)=response(comps{icomp,1},imode);            
        contourf(mesh.X,mesh.Y,real(Uplot),64,'linecolor','none');
%             countourfMultiDomain(mesh,real(Uplot));
        hold on;
        title(sprintf('Response %s, w=%.2f', comps{icomp,2},w))
        colorbar;
        daspect([1,1,1])

    subplot(5,2,2+(icomp-1)*2)
        Vplot = mesh.X*nan;
        Vplot(mesh.usedInd)=force(comps{icomp,1},imode);
        contourf(mesh.X,mesh.Y,real(Vplot),64,'linecolor','none');
%             countourfMultiDomain(mesh,real(Vplot));
        title(sprintf('Force %s, w=%.2f', comps{icomp,2},w))
        colorbar;
        daspect([1,1,1])
    end
end


