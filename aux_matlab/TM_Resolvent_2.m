function [S,V,fList,SS_conv] = TM_Resolvent_2(TM_setup,TM_setup_adj,deltaF,nf,n_Iter,tol,W,invW,B,C,mRSVD)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    if ~exist('mRSVD','var')     ; mRSVD=1; end
    
    dt        = TM_setup.dt     ;
    nq        = TM_setup.n      ; 
    nq_multi  = TM_setup.n_multi;
    
    %get block size
    if nf==1
        T_block     = dt     ;
        n_block     = round(T_block/dt);
        i_sampling  = n_block/nf;
        i_out       = 0;
        dt_sampling = 1;
    else
        T_block     = 1./abs(deltaF);
        n_block     = round(T_block/dt);

        D=find(mod(n_block, 1:n_block)==0); 
        [~,i]=min((D-nf)+1e99*(D<nf));
        nf=D(i);
        
        i_sampling  = n_block/nf;
        dt_sampling = dt*i_sampling;
        i_out       = (0:nf-1)*i_sampling;
    end
    
    if mod(nf,2)==0; fList = (-((nf/2))  :(nf/2)-1 )*deltaF;
    else           ; fList = (- (nf-1)/2 :(nf-1)/2)*deltaF;
    end
    fList=ifftshift(fList);
    
    % Alocate input and output vectors
    f_in          = zeros(nq,mRSVD,nf,n_Iter);
    f_out         = zeros(nq,mRSVD,nf,n_Iter);

    
    
    % Create initial random seeds
    f_prev(:,:,:) = randn(nq,mRSVD,nf)+randn(nq,mRSVD,nf)*1i;
    % Make seeds at each frequency ortogonal
    for i_f = 1:nf
        [f_prev(:,:,i_f),~] = qr(f_prev(:,:,i_f),0);
    end
    
    %Perform nIter iterations
    for i_iter = 1:n_Iter
        disp(['Computing iteration : ' num2str(i_iter)])
        %save input data
        f_in(:,:,:,i_iter) =  f_prev ; 

        %define forcing term for the direct problem
        for i_m=1:mRSVD                
            %%
            q0   = zeros(nq*nq_multi,1);
            
            %Apply the Resolvent operator
            flagAdj = false;
%             q_hat2 = CRB_TM(TM_setup,fList,q0,squeeze(f_prev(:,i_m,:)),B,C,i_out,n_block,inf,tol,flagAdj);
            
            q_hat = CRB_TM_gmres(TM_setup,fList,squeeze(f_prev(:,i_m,:)),B,C,i_out,n_block,tol,flagAdj);

%             norm(q_hat-q_hat2)
            q_hat = W*q_hat;
            
            
            % Normalize frequency components
            norms = zeros(nf,1);
            for i=1:nf
                norms(i)    = norm(q_hat(:,i))   ;
                q_hat(:,i)  = q_hat(:,i)/norms(i);
            end

            %Apply the Adjoint Resolvent operator
            flagAdj = true;
%             f_hat2 = CRB_TM      (TM_setup_adj,fList,q0,q_hat,B,C,i_out,n_block,inf,tol,flagAdj);
            f_hat = CRB_TM_gmres(TM_setup_adj,fList,   q_hat,B,C,i_out,n_block,tol,flagAdj);

            f_hat = invW*f_hat;
            
            % De-normalize frequenci components
            for i=1:nf
                f_hat(:,i)  = f_hat(:,i)*norms(i);
            end

            f_out(:,i_m,:,i_iter) = f_hat;

        end
        %% Make next inputs ortogonal to all previous ones
        for i_f = 1:nf
            f  = reshape( f_out(:,:,i_f,i_iter  ),nq,mRSVD       ) ;  % current outputs
            F  = reshape( f_in (:,:,i_f,1:i_iter),nq,mRSVD*i_iter) ;    %previous input
            f = f - F*(F'*f);  % keep only f that is orthogonal to F 
            [f,~]= qr(f,0);
            f_prev(:,:,i_f) = f;
        end
    end
    
    %% Compute gains based on the Arnold decomposition constructed
    S       = zeros(n_Iter*mRSVD,nf)      ; %gains 
    SS_conv = nan(n_Iter*mRSVD,n_Iter,nf) ; %gains convergence history

    for i_f = 1:nf
        F_in   = reshape( f_in  (:,:,i_f,:),nq,mRSVD*n_Iter) ;    % Inputs
        F_out  = reshape( f_out (:,:,i_f,:),nq,mRSVD*n_Iter) ;    % Output
        H = F_in'*F_out;
        
        [psi,SS]        = eig(H);
        SS              = sqrt(diag(SS));
        [S(:,i_f),order]= sort(SS,'descend');
        V(:,:,i_f)      = F_in*psi(:,order); 
        
        for j=1:n_Iter
            SS_conv(1:j*mRSVD,j,i_f) = sort(sqrt( eig(H(1:j*mRSVD,1:j*mRSVD))) ,'descend');
        end
    end
    
end

