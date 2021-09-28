function [gain,U_out,V_in] = ResAnalysis(L0,F,invF,omega,noEigs,B,C,FILTER_mesh,FILTER_ct_mesh)%,filter_f,r_1D,z_1D,X,Y,flagLowMemory)
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
    nDOFs=size(L0,1);
    nMesh=nDOFs/5;
    
    if ~exist('F'); f=speye(nDOFs);end
    if ~exist('invF'); invf=speye(nDOFs);end
    if ~exist('B'); B=speye(nDOFs);end
    if ~exist('C'); B=speye(nDOFs);end
    if ~exist('flagLowMemory'); flagLowMemory=true;end
    
    if ~exist('FILTER_mesh') & ~exist('FILTER_ct_mesh') 
        %no filter specified
        disp('No filter specified.')
        FILTER = @(x) x;
        FILTER_ct = @(x) x;
    elseif exist('FILTER_mesh')  + exist('FILTER_ct_mesh') == 1
        error('Only one of FILTER and FILTER_ct were ser. Both are needed.')
    else
        disp('Using speficified filters.')
        FILTER    = @(x) [  FILTER_mesh(x((1:nMesh)+nMesh*0))  ; 
                            FILTER_mesh(x((1:nMesh)+nMesh*1))  ;
                            FILTER_mesh(x((1:nMesh)+nMesh*2))  ;
                            FILTER_mesh(x((1:nMesh)+nMesh*3))  ;
                            FILTER_mesh(x((1:nMesh)+nMesh*4))   ] ;
                        
        FILTER_ct = @(x) [  FILTER_ct_mesh(x( (1:nMesh)+nMesh*0)) ; 
                            FILTER_ct_mesh(x( (1:nMesh)+nMesh*1)) ;
                            FILTER_ct_mesh(x( (1:nMesh)+nMesh*2)) ;
                            FILTER_ct_mesh(x( (1:nMesh)+nMesh*3)) ;
                            FILTER_ct_mesh(x( (1:nMesh)+nMesh*4))   ] ;
    end
      

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Input-output analysis %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % see Nichols, CTR report 2014 "Input-output analysis of high-speed jet noise"
    disp(['    calling intput-output analysis with omega=' num2str(omega)]);

    if size(F,1)==1 || size(F,2)==1
        F       = diag(sparse(F));
    end
    if size(invF,1)==1 || size(invF,2)==1
        invF    = diag(sparse(invF));
    end

    % get input and output matrices from user-defined routine
    
    % LU-decomposition of L-sigma*I
    tic
    OMEGA = omega*speye(nDOFs)%G+1i*spdiags(repmat(sponge,5,1),0,nDOFs,nDOFs);
%     LsI                 = L0-1i*omega*speye(nDOFs);
    LsI                 = L0-1i*OMEGA;
    [LL,UU,pp,qq,rr]    = lu(LsI);
    time = toc;
    disp(['    elapsed time - LU-decomposition of L-sigma*I: ' datestr(time/24/3600, 'HH:MM:SS')]);

   

    %%


%         HtrH    = @(v) invF*(B_ct*(rr_ct\(pp_ct*(LL_ct\(UU_ct\(qq_ct*(C_ct*(F*(C*(qq*(UU\(LL\(pp*(rr\(B*v)))))))))))))));
    H      = @(v) FILTER(C*(qq*(UU\(LL\(pp*(rr\(B*FILTER(v))))))));
    if flagLowMemory % saves memomy but slows down HTr
        LL_ct   = LL';    
        UU_ct   = UU';    
        pp_ct   = pp';    qq_ct   = qq'; rr_ct   = rr';
        C_ct    = C';    B_ct    = B';
        Htr    = @(v) FILTER_ct(B_ct*(rr_ct\(pp_ct*(LL_ct\(UU_ct\(qq_ct*(C_ct*FILTER_ct(v))))))));
    else % higher memory footprint, but faster interations of HTr
        Htr    = @(v) FILTER_ct(B'*(rr'\(pp'*(LL'\(UU'\(qq'*(C'*FILTER_CT(v))))))));
    end

    %% adapts function for running on Octave or Matlab
    isOctave =  exist('OCTAVE_VERSION', 'builtin') ~= 0;
    if isOctave
        %Octave does not hangle complex valued functions.
        comp2real = @(x) [real(x);imag(x)];
        real2comp = @(x) real(x(1:end/2))+1i*x(end/2+1:end);
        
        HtrH   = @(v) comp2real(invF*Htr(F*H(real2comp(v))));
    else    
        HtrH   = @(v) invF*Htr(F*H(v));
    end
    
    opts.tol            = eps;
    opts.disp           = 3;
    opts.issym = false;

    % input
    tic
    if isOctave
        [V_in_tmp, OMEGA]       = eigs(HtrH, 2*nDOFs, noEigs, 'lm', opts);
        for i=1:size(V_in_tmp,2)
            V_in(:,i) = real2comp(V_in_tmp(:,i));
        end
    else
         [V_in, OMEGA]       = eigs(HtrH, nDOFs, noEigs, 'lm', opts);
    end
    gain                = sqrt(real(diag(OMEGA))); % singular values of H given by sqrt. of eigenvalues of H*H, should be real (imag. part is O(eps))

    time = toc;
    disp(['    elapsed time - SVD via ''eigs'' using matrix-free Arnoldi: ' datestr(time/24/3600, 'HH:MM:SS')]);            %%

    % output
    U_out = zeros(size(V_in));
    for eig_i = 1:noEigs
        % normalization in inner product norm
        V_in(:,eig_i)  = V_in(:,eig_i)/sqrt(V_in(:,eig_i)'*F*V_in(:,eig_i));
        U_out(:,eig_i) = H(V_in(:,eig_i))/gain(eig_i);
    end
