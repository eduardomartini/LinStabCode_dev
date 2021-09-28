function [DR,D2R,DZ,D2Z,D2RZ,D2ZR] = CreateDiffMatrices_Axy(mesh,m)
    % [DR,D2R,DZ,D2Z] = CreateDiffMatrices_Axy(mesh,m)
    % Create global differentiatial matricies to be used in an axysimmetric
    % analysis. 
tic
pipeBC=false;

pipe_idx_bottom     = [];
pipe_idx_top        = [];

% [Nr,Nz] = size(mesh.X);
ndofs=size(mesh.DW.Dx,1);

z = mesh.X;
r = mesh.Y;

Dz = mesh.DW.Dx;
Dr = mesh.DW.Dy;
Dr_symm  = mesh.DW.Dy_symm;       
Dr_asymm = mesh.DW.Dy_asymm;     

D2z         = mesh.DW.D2x;
D2r         = mesh.DW.D2y;
D2r_symm    = mesh.DW.D2y_symm;       
D2r_asymm   = mesh.DW.D2y_asymm;     

D2zr         = mesh.DW.Dxy;
D2rz         = mesh.DW.Dyx;
D2zr_symm    = mesh.DW.Dxy_symm;       
D2rz_symm    = mesh.DW.Dyx_symm;       
D2rz_asymm   = mesh.DW.Dxy_asymm;       
D2zr_asymm   = mesh.DW.Dyx_asymm;       


%%
calcMethod='';
Z = 0*speye(ndofs,ndofs); % need a zero diagonal sparse matrix... is there a better way?

% r dependent derivatives (depend on symmetry)
if  m==0
    DR   = blkdiag(Dr_symm,Dr_asymm,Dr_asymm,Dr_symm,Dr_symm);
    D2R  = blkdiag(D2r_symm,D2r_asymm,D2r_asymm,D2r_symm,D2r_symm);
    D2RZ = blkdiag(D2rz_symm,D2rz_asymm,D2rz_asymm,D2rz_symm,D2rz_symm);
    D2ZR = blkdiag(D2zr_symm,D2zr_asymm,D2zr_asymm,D2zr_symm,D2zr_symm);
elseif m==1
    DR   = blkdiag(Dr_asymm,Dr_symm,Dr_symm,Dr_asymm,Dr_asymm);
    D2R  = blkdiag(D2r_asymm,D2r_symm,D2r_symm,D2r_asymm,D2r_asymm);
    D2RZ = blkdiag(D2rz_asymm,D2rz_symm,D2rz_symm,D2rz_asymm,D2rz_asymm);
    D2ZR = blkdiag(D2zr_asymm,D2zr_symm,D2zr_symm,D2zr_asymm,D2zr_asymm);
elseif m>=2
    DR   = blkdiag(Dr_asymm,Dr_asymm,Dr_asymm,Dr_asymm,Dr_asymm);
    D2R  = blkdiag(D2r_asymm,D2r_asymm,D2r_asymm,D2r_asymm,D2r_asymm);
    D2RZ = blkdiag(D2rz_asymm,D2rz_asymm,D2rz_asymm,D2rz_asymm,D2rz_asymm);
    D2ZR = blkdiag(D2zr_asymm,D2zr_asymm,D2zr_asymm,D2zr_asymm,D2zr_asymm);
end

% z dependent derivatives (do not depend on symmetry)
DZ      = kron(speye(5,5),Dz);
D2Z     = kron(speye(5,5),D2z);

clear Z

disp(['    elapsed time - Differentiation matrices: ' datestr(toc/24/3600, 'HH:MM:SS')]);
