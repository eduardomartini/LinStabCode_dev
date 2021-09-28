function [defMesh] = DeformMesh(mesh,X2,Y2)
    % newMesh = DeformMesh(mesh,X2,Y2)
    % Deforms the 'mesh' object to the coordinates given by X2 and Y2.
    % Derivative matrices and integration weights are updated accordinly. 
    if isfield(mesh,'originalMesh')
        mesh=mesh.originalMesh;
    end
    
    defMesh                 = mesh;
    defMesh.X               = X2;
    defMesh.Y               = Y2;
    defMesh.originalMesh    = mesh;
    
    X1 = mesh.X;
    Y1 = mesh.Y;
    
    
    
    mDx  = mesh.DW.Dx;
    mD2x = mesh.DW.D2x;
    mDy  = mesh.DW.Dy;
    mD2y = mesh.DW.D2y;
    
    %% Construct new mesh and prepare diff matrices to be transformed

    ngp = mesh.ngp;
    
    X2 = X2(mesh.usedInd);
    Y2 = Y2(mesh.usedInd);
    
    
    Dx          = mesh.DW.Dx;
    D2x         = mesh.DW.D2x;
    Dy          = mesh.DW.Dy;
    Dy_symm     = mesh.DW.Dy_symm;
    Dy_asymm    = mesh.DW.Dy_asymm;
    D2y         = mesh.DW.D2y;
    D2y_symm    = mesh.DW.D2y_symm;
    D2y_asymm   = mesh.DW.D2y_asymm;
    Dxy         = mesh.DW.Dxy;
    Dxy_symm    = mesh.DW.Dxy_symm;
    Dxy_asymm   = mesh.DW.Dxy_asymm;
    Dyx         = mesh.DW.Dyx;
    Dyx_symm    = mesh.DW.Dyx_symm;
    Dyx_asymm   = mesh.DW.Dyx_asymm;
        
    %% Compute transformation derivatives
    %     Jacobian
    J = zeros(2,2,ngp);
    J(1,1,:) = mDx*X2(:);
    J(1,2,:) = mDx*Y2(:);

    J(2,1,:) = mDy*X2(:);
    J(2,2,:) = mDy*Y2(:);             
    % 
    dJ = zeros(2,2,ngp,2);
    dJ(1,1,:,1) = mD2x*X2(:);
    dJ(1,2,:,1) = mD2x*Y2(:);

    dJ(2,1,:,1) = mDx*(mDy*X2(:));
    dJ(2,2,:,1) = mDx*(mDy*Y2(:)); 

    dJ(1,1,:,2) = mDy*(mDx*X2(:));
    dJ(1,2,:,2) = mDy*(mDx*Y2(:));

    dJ(2,1,:,2) = mD2y*X2(:);
    dJ(2,2,:,2) = mD2y*Y2(:); 

    % Inverse transform derivative
    Ji  = zeros(2,2,ngp);
    det = J(1,1,:).*J(2,2,:)- J(2,1,:).*J(1,2,:);
    Ji(1,1,:)  =  J(2,2,:)./det;
    Ji(2,1,:)  = -J(2,1,:)./det;
    Ji(1,2,:)  = -J(1,2,:)./det;
    Ji(2,2,:)  =  J(1,1,:)./det;
        

    

    %% Compute new derivative matrices
    %% First Derivatives
    % f(xi)=f(xi(sigj))
    % df/dxi = df/dsigj dsigj/dxi  =  df/dsigj (dxj/dsigi)^-1  = Ji_ij * df/dsigj

    Dx_new          = spdiags(squeeze(Ji(1,1,:)),0,ngp,ngp)*Dx ...
                    + spdiags(squeeze(Ji(1,2,:)),0,ngp,ngp)*Dy;
    Dy_new          = spdiags(squeeze(Ji(2,1,:)),0,ngp,ngp)*Dx ...
                    + spdiags(squeeze(Ji(2,2,:)),0,ngp,ngp)*Dy;
    Dy_symm_new     = spdiags(squeeze(Ji(2,1,:)),0,ngp,ngp)*Dx ...
                    + spdiags(squeeze(Ji(2,2,:)),0,ngp,ngp)*Dy_symm;
    Dy_asymm_new    = spdiags(squeeze(Ji(2,1,:)),0,ngp,ngp)*Dx ...
                    + spdiags(squeeze(Ji(2,2,:)),0,ngp,ngp)*Dy_asymm;

    %% Second Derivatives
    % ddf/dxidxj = Ji_ik  ddf/dsigkdsigl Ji_jl  +         d(Ji_ik)/dxj   * df/dsigk 
    %            = Ji_ik  ddf/dsigkdsigl Ji_jl  + Ji_jl * d(Ji_ik)/dsigl * df/dsigk 
    %            = Ji_ik  ddf/dsigkdsigl Ji_jl  + Q_ijk                  * df/dsigk 

    JiT = Ji;
    JiT(1,2,:) = Ji(2,1,:);
    JiT(2,1,:) = Ji(1,2,:);
    for i=1:4
        A=zeros(2,2); A(i)=1;
        a(:,:,:,i)  =  multiprod(JiT,multiprod(A, Ji));
    end
    axx=squeeze(a(1,1,:,:));
    axy=squeeze(a(1,2,:,:));
    ayx=squeeze(a(2,1,:,:));
    ayy=squeeze(a(2,2,:,:));

    % d(Ji)/dsigl = - Ji*dJ/dsigl*Ji 

    dJi(:,:,:,1)=-multiprod(multiprod(Ji,dJ(:,:,:,1)),Ji);
    dJi(:,:,:,2)=-multiprod(multiprod(Ji,dJ(:,:,:,2)),Ji);

    % Q_ijk = Ji_jl * d(Ji_ik)/dsigl
    Q = zeros(2,2,2,ngp);
    for i=1:2
        for j=1:2
            for k=1:2
                Q(i,j,k,:) = Ji(j,1,:).*dJi(i,k,:,1) + Ji(j,2,:).*dJi(i,k,:,2);
            end
        end    
    end

    Dxy = Dx*Dy;
    Dyx = Dy*Dx;

    diags = @(x)spdiags(x,0,ngp,ngp);

     D2x_new =  diags(axx(:,1))*D2x + diags(axy(:,1))*Dxy + ...
                diags(ayx(:,1))*Dyx + diags(ayy(:,1))*D2y + ...
                diags(squeeze(Q(1,1,1,:)))*Dx + diags(squeeze(Q(1,1,2,:)))*Dy ;
     Dxy_new =  diags(axx(:,2))*D2x + diags(axy(:,2))*Dxy + ...
                diags(ayx(:,2))*Dyx + diags(ayy(:,2))*D2y + ...
                diags(squeeze(Q(1,2,1,:)))*Dx + diags(squeeze(Q(1,2,2,:)))*Dy ;
     Dyx_new =  diags(axx(:,3))*D2x + diags(axy(:,3))*Dxy + ...
                diags(ayx(:,3))*Dyx + diags(ayy(:,3))*D2y + ...
                diags(squeeze(Q(2,1,1,:)))*Dx + diags(squeeze(Q(2,1,2,:)))*Dy ;
     D2y_new =  diags(axx(:,4))*D2x + diags(axy(:,4))*Dxy + ...
                diags(ayx(:,4))*Dyx + diags(ayy(:,4))*D2y + ...
                diags(squeeze(Q(2,2,1,:)))*Dx + diags(squeeze(Q(2,2,2,:)))*Dy ;
     
     Dxy_symm_new =  diags(axx(:,2))*D2x            + diags(axy(:,2))*Dxy_symm + ...
                     diags(ayx(:,2))*Dyx_symm       + diags(ayy(:,2))*D2y_symm + ...
                     diags(squeeze(Q(1,2,1,:)))*Dx  + diags(squeeze(Q(1,2,2,:)))*Dy_symm ;
     Dyx_symm_new =  diags(axx(:,3))*D2x            + diags(axy(:,3))*Dxy_symm + ...
                     diags(ayx(:,3))*Dyx_symm       + diags(ayy(:,3))*D2y_symm + ...
                     diags(squeeze(Q(2,1,1,:)))*Dx  + diags(squeeze(Q(2,1,2,:)))*Dy_symm ;
     D2y_symm_new =  diags(axx(:,4))*D2x            + diags(axy(:,4))*Dxy_symm + ...
                     diags(ayx(:,4))*Dyx_symm       + diags(ayy(:,4))*D2y_symm + ...
                     diags(squeeze(Q(2,2,1,:)))*Dx  + diags(squeeze(Q(2,2,2,:)))*Dy_symm ;
     
     Dxy_asymm_new=  diags(axx(:,2))*D2x            + diags(axy(:,2))*Dxy_asymm + ...
                     diags(ayx(:,2))*Dyx_asymm       + diags(ayy(:,2))*D2y_asymm + ...
                     diags(squeeze(Q(1,2,1,:)))*Dx  + diags(squeeze(Q(1,2,2,:)))*Dy_asymm ;
     Dyx_asymm_new=  diags(axx(:,3))*D2x            + diags(axy(:,3))*Dxy_asymm + ...
                     diags(ayx(:,3))*Dyx_asymm       + diags(ayy(:,3))*D2y_asymm + ...
                     diags(squeeze(Q(2,1,1,:)))*Dx  + diags(squeeze(Q(2,1,2,:)))*Dy_asymm ;
     D2y_asymm_new=  diags(axx(:,4))*D2x            + diags(axy(:,4))*Dxy_asymm + ...
                     diags(ayx(:,4))*Dyx_asymm       + diags(ayy(:,4))*D2y_asymm + ...
                     diags(squeeze(Q(2,2,1,:)))*Dx  + diags(squeeze(Q(2,2,2,:)))*Dy_asymm ;
                  
            
    % update intregration weights
    Wnew        = zeros(size(mesh.DW.W));
    Wnew_symm   = zeros(size(mesh.DW.W));
    Wnew_asymm  = zeros(size(mesh.DW.W));
    Wnew(mesh.usedInd)       = mesh.DW.W(mesh.usedInd).*abs(det(:));
    Wnew_symm(mesh.usedInd)  = mesh.DW.W_symm(mesh.usedInd).*abs(det(:));
    Wnew_asymm(mesh.usedInd) = mesh.DW.W_asymm(mesh.usedInd).*abs(det(:));

            
    defDW.Dx        = Dx_new ;
    defDW.Dy        = Dy_new ;
    defDW.Dy_symm   = Dy_symm_new ;
    defDW.Dy_asymm  = Dy_asymm_new ;

    defDW.D2x       = D2x_new ;
    defDW.Dxy       = Dxy_new ;
    defDW.Dyx       = Dyx_new ;
    defDW.D2y       = D2y_new ;
    
    defDW.Dxy_symm  = Dxy_symm_new ;
    defDW.Dyx_symm  = Dyx_symm_new ;
    defDW.D2y_symm  = D2y_symm_new ;

    defDW.Dxy_asymm = Dxy_asymm_new ;
    defDW.Dyx_asymm = Dyx_asymm_new ;
    defDW.D2y_asymm = D2y_asymm_new ;

    defDW.W         = Wnew;
    defDW.W_symm    = Wnew_symm;
    defDW.W_asymm   = Wnew_asymm;
    
    defMesh.DW              = defDW;

