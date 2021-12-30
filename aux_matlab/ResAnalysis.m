function [gain,U_out,V_in] = ResAnalysis(LinProb,Wen,invWen,omega,noEigs,B,C)
    % [gain,U_out,V_in] = ResAnalysis(LHS,RHS,F,invF,omega,noEigs,B,C)%,filter_f,r_1D,z_1D,X,Y)
    % Performs classical resolvent analysis using "eigs" function of the
    % linear operator R=(L0-i*omega)^-1 using input and output matrices 
    % given by B and C. The adjoint is obtained with respect to the norm F, 
    % and Finv = F^-1. 
    % Inputs : 
    %     noEigs : number of eigen values desired.
    % Outputs :
    %     gain  : Resolvent gains
    %     U_out : optimal outputs (Responses)
    %     V_in  : optimal inputs (Forcing)

    % reduce to generalized to standard EVP

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Input-output analysis %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % see Nichols, CTR report 2014 "Input-output analysis of high-speed jet noise"
    disp(['    calling intput-output analysis with omega=' num2str(omega)]);
    
    % LU-decomposition of L-sigma*I
    tic
    
    A    = LinProb.A;
    T_exp = LinProb.T_exp;
    T_cont= LinProb.T_cont;
    
    A_red   = T_cont*A*T_exp;
    
    I       = speye(size(A_red));
    LsI     = A_red-1i*omega*I;
    
    [LL,UU,pp,qq,rr]    = lu(LsI);
    
    time = toc;
    disp(['    elapsed time - LU-decomposition of L-sigma*I: ' datestr(time/24/3600, 'HH:MM:SS')]);

   

    %%
        if size(B) == size(A)
        B       = T_cont*B*T_exp;
    end
    if size(C) == size(A)
        C       = T_cont*C*T_exp;
    end
    if size(Wen) == size(A)
        Wen     = T_cont*Wen*T_exp;
    end
    if size(invWen) == size(A)
        invWen  = T_cont*invWen*T_exp;
    end
    
    LL_ct   = LL';    
    UU_ct   = UU';    
    pp_ct   = pp';    
    qq_ct   = qq'; 
    rr_ct   = rr';
    C_ct    = C';    
    B_ct    = B';
        
    FILTER    = @(x) T_cont*LinProb.FILTER   (T_exp  *x);
    FILTER_ct = @(x) T_exp'*LinProb.FILTER_ct(T_cont'*x);
    
    %Adjust size of B and C to match the problem with imposed B.C.

    
    
    R      = @(v) FILTER   (qq   *(UU   \(LL   \(pp   *(   rr\FILTER   (v))))));
    R_ct   = @(v) FILTER_ct(rr_ct\(pp_ct*(LL_ct\(UU_ct\(qq_ct*FILTER_ct(v))))));

    
    CRB    = @(v) C   *R   (B   *v);
    CRB_ct = @(v) B_ct*R_ct(C_ct*v);
    
    %% adapts function for running on Octave or Matlab
      
    RadjR   = @(v) invWen*CRB_ct(Wen*CRB(v));
    
    opts.tol            = eps;
    opts.disp           = 3;
    opts.issym = false;

    % input
    tic

    [V_in, GAIN2]       = eigs(RadjR, size(A_red,1), noEigs, 'lm', opts);
    
    gain                = sqrt(real(diag(GAIN2))); % singular values of H given by sqrt. of eigenvalues of H*H, should be real (imag. part is O(eps))

    time = toc;
    disp(['    elapsed time - SVD via ''eigs'' using matrix-free Arnoldi: ' datestr(time/24/3600, 'HH:MM:SS')]);            %%

    % output
    U_out = zeros(size(V_in));
    for eig_i = 1:noEigs
        % normalization in inner product norm
        V_in(:,eig_i)  = V_in(:,eig_i)/sqrt(V_in(:,eig_i)'*Wen*V_in(:,eig_i));
        U_out(:,eig_i) = R(V_in(:,eig_i))/gain(eig_i);
    end
    
    V_in  = T_exp*V_in;
    U_out = T_exp*U_out;
