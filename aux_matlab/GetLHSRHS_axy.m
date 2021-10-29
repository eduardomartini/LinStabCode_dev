function [LHS,RHS,indexes] = GetLHSRHS_axy(mesh,BF,kzm,Re,sponge,model)
    % [LHS,RHS] = GetLHSRHS(mesh,BF,m,Re,sponge)
    % Creates LHS and RHS operators of the LNS equations, such that
    % LHS dq/dt = RHS q , q=[rho,u,v,w,T].
    % Inputs ; 
    %     mesh : mesh object containing the spatial coordiantes,
    %     diferentiation matrices and integration weights.
    %     BF   : object containing the baseflow to be used
    %     sponge : sponge function to be used.
    % Outputs :
    %     LHS and RHS matrices.
    
    tic
    % Get mesh and scalar diff. matrices.
    z = mesh.X(mesh.usedInd);
    Dz = mesh.DW.Dx;
    D2z = mesh.DW.D2x;
    
    r = mesh.Y(mesh.usedInd);
    Dr = mesh.DW.Dy;
    D2r = mesh.DW.D2y;

    Drz = mesh.DW.Dxy;
    Dzr = mesh.DW.Dyx;

        
    nGridPoints = mesh.ngp;
    
    % Get meanflow in vector form
    U   = BF.U(mesh.usedInd);
    V   = BF.V(mesh.usedInd);
    W   = BF.W(mesh.usedInd);
    RHO = BF.RHO(mesh.usedInd);
    T   = BF.T(mesh.usedInd);
    MU  = BF.MU(mesh.usedInd);
    P   = RHO.*T./(BF.Ma^2*BF.kappa);
    dmudT    = BF.dmudT(mesh.usedInd);
    d2mudT2  = BF.d2mudT2(mesh.usedInd);
    dMUdT    = BF.dMUdT(mesh.usedInd);
    d2MUdT2  = BF.d2MUdT2(mesh.usedInd);

    cv   = BF.cv;
    c1   = BF.c1;
    c2   = BF.c2;
    kappa= BF.kappa;

    %aux ùatrices
    Z = zeros(nGridPoints, 1);
    I = ones(nGridPoints, 1);
    R = reshape(r, nGridPoints, 1);

    %Sponge matrix
    Asponge = spdiags( [sponge(:);sponge(:);sponge(:);sponge(:);sponge(:)],0,nGridPoints*5,nGridPoints*5);
        
    %% COMPUTE BASEFLOW DERIVATIVES
    tic
    % Base flow derivatives
    dUdr    = Dr*U;     d2Udr2    = D2r*U;
    dVdr    = Dr*V;     d2Vdr2    = D2r*V;
    dWdr    = Dr*W;     d2Wdr2    = D2r*W;
    dRHOdr  = Dr*RHO;   d2RHOdr2  = D2r*RHO;
    dTdr    = Dr*T;     d2Tdr2    = D2r*T;
    dUdz    = Dz*U;     d2Udz2    = D2z*U;
    dVdz    = Dz*V;     d2Vdz2    = D2z*V;
    dWdz    = Dz*W;     d2Wdz2    = D2z*W;
    dRHOdz  = Dz*RHO;   d2RHOdz2  = D2z*RHO;
    dTdz    = Dz*T;     d2Tdz2    = D2z*T;
    d2Udrz  = Dz*dUdr;
    d2Vdrz  = Dz*dVdr;
    d2Wdrz  = Dz*dWdr; 
    dPdr    = Dr*P;
    dPdz    = Dz*P;
    
    
    %% Coefficient matrices
    %get global diff. matrices
%     D2ZR=DR*DZ;
    
    %get coefficients
    if strcmp(model,'2D')
        Ma=BF.Ma;
        Pr=BF.Pr;
        kz = kzm;
        NrNz=nGridPoints;
        DR  =blkdiag(Dr ,Dr ,Dr ,Dr ,Dr );
        D2R =blkdiag(D2r,D2r,D2r,D2r,D2r);
        DZ  =blkdiag(Dz ,Dz ,Dz ,Dz ,Dz );
        D2Z =blkdiag(D2z,D2z,D2z,D2z,D2z);
        D2RZ=blkdiag(Dzr,Dzr,Dzr,Dzr,Dzr);
        getCoeffsLaminar_Cartesian
    else
        m=kzm;
        [DR,D2R,DZ,D2Z,D2RZ,D2ZR] =   CreateDiffMatrices_Axy(mesh,m);
        getCoeffsLaminar
    end

    %build matrices from coefficients
    A0 =    ...
        [   ...
        diag(sparse(A0_11)) diag(sparse(A0_12)) diag(sparse(A0_13)) diag(sparse(A0_14)) diag(sparse(A0_15))
        diag(sparse(A0_21)) diag(sparse(A0_22)) diag(sparse(A0_23)) diag(sparse(A0_24)) diag(sparse(A0_25))
        diag(sparse(A0_31)) diag(sparse(A0_32)) diag(sparse(A0_33)) diag(sparse(A0_34)) diag(sparse(A0_35))
        diag(sparse(A0_41)) diag(sparse(A0_42)) diag(sparse(A0_43)) diag(sparse(A0_44)) diag(sparse(A0_45))
        diag(sparse(A0_51)) diag(sparse(A0_52)) diag(sparse(A0_53)) diag(sparse(A0_54)) diag(sparse(A0_55))
        ];

    Ar =    ...
        [   ...
        diag(sparse(Ar_11)) diag(sparse(Ar_12)) diag(sparse(Ar_13)) diag(sparse(Ar_14)) diag(sparse(Ar_15))
        diag(sparse(Ar_21)) diag(sparse(Ar_22)) diag(sparse(Ar_23)) diag(sparse(Ar_24)) diag(sparse(Ar_25))
        diag(sparse(Ar_31)) diag(sparse(Ar_32)) diag(sparse(Ar_33)) diag(sparse(Ar_34)) diag(sparse(Ar_35))
        diag(sparse(Ar_41)) diag(sparse(Ar_42)) diag(sparse(Ar_43)) diag(sparse(Ar_44)) diag(sparse(Ar_45))
        diag(sparse(Ar_51)) diag(sparse(Ar_52)) diag(sparse(Ar_53)) diag(sparse(Ar_54)) diag(sparse(Ar_55))
        ];

    Az =    ...
        [   ...
        diag(sparse(Az_11)) diag(sparse(Az_12)) diag(sparse(Az_13)) diag(sparse(Az_14)) diag(sparse(Az_15))
        diag(sparse(Az_21)) diag(sparse(Az_22)) diag(sparse(Az_23)) diag(sparse(Az_24)) diag(sparse(Az_25))
        diag(sparse(Az_31)) diag(sparse(Az_32)) diag(sparse(Az_33)) diag(sparse(Az_34)) diag(sparse(Az_35))
        diag(sparse(Az_41)) diag(sparse(Az_42)) diag(sparse(Az_43)) diag(sparse(Az_44)) diag(sparse(Az_45))
        diag(sparse(Az_51)) diag(sparse(Az_52)) diag(sparse(Az_53)) diag(sparse(Az_54)) diag(sparse(Az_55))
        ];

    Arz =   ...
        [   ...
        diag(sparse(Arz_11)) diag(sparse(Arz_12)) diag(sparse(Arz_13)) diag(sparse(Arz_14)) diag(sparse(Arz_15))
        diag(sparse(Arz_21)) diag(sparse(Arz_22)) diag(sparse(Arz_23)) diag(sparse(Arz_24)) diag(sparse(Arz_25))
        diag(sparse(Arz_31)) diag(sparse(Arz_32)) diag(sparse(Arz_33)) diag(sparse(Arz_34)) diag(sparse(Arz_35))
        diag(sparse(Arz_41)) diag(sparse(Arz_42)) diag(sparse(Arz_43)) diag(sparse(Arz_44)) diag(sparse(Arz_45))
        diag(sparse(Arz_51)) diag(sparse(Arz_52)) diag(sparse(Arz_53)) diag(sparse(Arz_54)) diag(sparse(Arz_55))
        ];

    Arr =   ...
        [   ...
        diag(sparse(Arr_11)) diag(sparse(Arr_12)) diag(sparse(Arr_13)) diag(sparse(Arr_14)) diag(sparse(Arr_15))
        diag(sparse(Arr_21)) diag(sparse(Arr_22)) diag(sparse(Arr_23)) diag(sparse(Arr_24)) diag(sparse(Arr_25))
        diag(sparse(Arr_31)) diag(sparse(Arr_32)) diag(sparse(Arr_33)) diag(sparse(Arr_34)) diag(sparse(Arr_35))
        diag(sparse(Arr_41)) diag(sparse(Arr_42)) diag(sparse(Arr_43)) diag(sparse(Arr_44)) diag(sparse(Arr_45))
        diag(sparse(Arr_51)) diag(sparse(Arr_52)) diag(sparse(Arr_53)) diag(sparse(Arr_54)) diag(sparse(Arr_55))
        ];

    Azz =   ...
        [   ...
        diag(sparse(Azz_11)) diag(sparse(Azz_12)) diag(sparse(Azz_13)) diag(sparse(Azz_14)) diag(sparse(Azz_15))
        diag(sparse(Azz_21)) diag(sparse(Azz_22)) diag(sparse(Azz_23)) diag(sparse(Azz_24)) diag(sparse(Azz_25))
        diag(sparse(Azz_31)) diag(sparse(Azz_32)) diag(sparse(Azz_33)) diag(sparse(Azz_34)) diag(sparse(Azz_35))
        diag(sparse(Azz_41)) diag(sparse(Azz_42)) diag(sparse(Azz_43)) diag(sparse(Azz_44)) diag(sparse(Azz_45))
        diag(sparse(Azz_51)) diag(sparse(Azz_52)) diag(sparse(Azz_53)) diag(sparse(Azz_54)) diag(sparse(Azz_55))
        ];

    ZZ = 0*speye(nGridPoints,nGridPoints); % need a zero diagonal sparse matrix... is there a better way?

    RHS =     ...
        [   ...
        diag(sparse(B_11)) diag(sparse(B_12)) diag(sparse(B_13)) diag(sparse(B_14)) diag(sparse(B_15))
        diag(sparse(B_21)) diag(sparse(B_22)) diag(sparse(B_23)) diag(sparse(B_24)) diag(sparse(B_25))
        diag(sparse(B_31)) diag(sparse(B_32)) diag(sparse(B_33)) diag(sparse(B_34)) diag(sparse(B_35))
        diag(sparse(B_41)) diag(sparse(B_42)) diag(sparse(B_43)) diag(sparse(B_44)) diag(sparse(B_45))
        diag(sparse(B_51)) diag(sparse(B_52)) diag(sparse(B_53)) diag(sparse(B_54)) diag(sparse(B_55))
        ];

    %% Construct the operator
%     LHS     = A0 + Ar*DR + Az*DZ + Arr*D2R + Azz*D2Z + Arz*DR*DZ - Asponge;
    LHS     = A0 + Ar*DR + Az*DZ + Arr*D2R + Azz*D2Z + Arz*D2RZ;% - Asponge*(RHS/1i);
    
    %% Construct the varibles and boundary indexes.
    % Row indices for primitive variables
    indexes = struct();

%     indexes.li = indexes.left;
%     indexes.ri = indexes.right; 
%     indexes.bi = indexes.bottom;
%     indexes.ti = indexes.top;

%     indexes.ti_rho  = indexes.ti;
%     indexes.bi_rho  = indexes.bi;
%     indexes.ri_rho  = indexes.ri;
%     indexes.li_rho  = indexes.li;
% 
%     indexes.ti_u    = indexes.ti   + nGridPoints;
%     indexes.bi_u    = indexes.bi   + nGridPoints;
%     indexes.ri_u    = indexes.ri   + nGridPoints;
%     indexes.li_u    = indexes.li   + nGridPoints;
% 
%     indexes.ti_v    = indexes.ti   + 2*nGridPoints;
%     indexes.bi_v    = indexes.bi   + 2*nGridPoints;
%     indexes.ri_v    = indexes.ri   + 2*nGridPoints;
%     indexes.li_v    = indexes.li   + 2*nGridPoints;
% 
%     indexes.ti_w    = indexes.ti   + 3*nGridPoints;
%     indexes.bi_w    = indexes.bi   + 3*nGridPoints;
%     indexes.ri_w    = indexes.ri   + 3*nGridPoints;
%     indexes.li_w    = indexes.li   + 3*nGridPoints;
% 
%     indexes.ti_T    = indexes.ti   + 4*nGridPoints;
%     indexes.bi_T    = indexes.bi   + 4*nGridPoints;
%     indexes.ri_T    = indexes.ri   + 4*nGridPoints;
%     indexes.li_T    = indexes.li   + 4*nGridPoints;
% 
%     indexes.ti_all  = [indexes.ti_rho,indexes.ti_T, ...
%                        indexes.ti_u,indexes.ti_v,indexes.ti_w];
%     indexes.bi_all  = [indexes.bi_rho,indexes.bi_T, ...
%                        indexes.bi_u,indexes.bi_v,indexes.bi_w];
%     indexes.ri_all  = [indexes.ri_rho,indexes.ri_T, ...
%                        indexes.ri_u,indexes.ri_v,indexes.ri_w];
%     indexes.li_all  = [indexes.li_rho,indexes.li_T, ...
%                        indexes.li_u,indexes.li_v,indexes.li_w];

    
    % Column indices for conserved variables
    indexes.rho_j   =        1        :   nGridPoints;
    indexes.u_j     =   nGridPoints+1 : 2*nGridPoints;
    indexes.v_j     = 2*nGridPoints+1 : 3*nGridPoints;
    indexes.w_j     = 3*nGridPoints+1 : 4*nGridPoints;
    indexes.T_j     = 4*nGridPoints+1 : 5*nGridPoints;
    
% 
%     time    = toc;
%     disp(['    elapsed time - Coefficient matrices: ' datestr(time/24/3600, 'HH:MM:SS')]);
