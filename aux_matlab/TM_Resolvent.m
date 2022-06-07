function [S,V,fList,SS_conv] = TM_Resolvent(TM_setup,TM_setup_adj,deltaF,nf,n_Iter,tol,W,invW,B,C)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    m         = 1               ;
    dt        = TM_setup.dt     ;
    nq        = TM_setup.n      ; 
    nq_multi  = TM_setup.n_multi;
    
%     nf  = length(fList);
    if mod(nf,2)==0; fList = (-((nf/2))  :(nf/2)-1 )*deltaF;
    else           ; fList = (- (nf-1)/2 :(nf-1)/2)*deltaF;
    end
    fList=ifftshift(fList);
    %get block size
    if nf==1
        T_block = dt     ;
        n_block     = round(T_block/dt);
        i_sampling  = n_block/nf;
        i_out       = 0;
        dt_sampling = 1;
    else
        T_block     = fList(2)-fList(1);
        n_block     = round(T_block/dt);
        i_sampling  = n_block/nf;
        dt_sampling = dt*i_sampling;
        i_out       = (0:nf-1)*i_sampling;
    end
    
    % Alocate input and output vectors
    f_in          = zeros(nq,m,nf,n_Iter);
    f_out         = zeros(nq,m,nf,n_Iter);

        
    
    % Create initial random seeds
    f_prev(:,:,:) = randn(nq,m,nf)+randn(nq,m,nf)*1i;
    % Make seeds at each frequency ortogonal
    for i_f = 1:nf
        [f_prev(:,:,i_f),~] = qr(f_prev(:,:,i_f),0);
    end
    
    %Perform nIter iterations
    for i_iter = 1:n_Iter
        %save input data
        f_in(:,:,:,i_iter) =  f_prev ; 

        %define forcing term for the direct problem
        for i_m=1:m                
            %%
            f = @(t) B*squeeze(f_prev(:,i_m,:))*exp(-2i*pi*fList.'*t);
            q0   = zeros(nq*nq_multi,1);
            
            [q,qhist]    = TM(q0,f ,inf,TM_setup,i_out,n_block,tol);
            %%
            q     = C*q(1:nq,:); 
            q_hat = ifft(q,[],2)*nf*dt_sampling;            
            q_hat = W*q_hat;
            
            
%             subplot(211)
% %                 [q_hat,-L0\squeeze(f_prev(:,i_m,:))]
% % 
%                 fadj = @(t) q_hat             *exp(-2i*pi*fList.'*t);
%                 tt = (0:size(qhist,2)-1)*dt;
%                 qtest = fadj(tt );
%                 plot(tt,qhist(1,:),tt,qtest(1,:),':',[i_out+n_block*1]*dt,q(1,:),'o');

            %define forcing term for the adjoint problem
            
            % Normalize frequenci components
            norms = zeros(nf,1);
            for i=1:nf
                norms(i)    = norm(q_hat(:,i))   ;
                q_hat(:,i)  = q_hat(:,i)/norms(i);
            end

            
            fadj = @(t) C*q_hat             *exp(2i*pi*fList.'*t);
            [f,fhist]   = TM(q0,fadj ,inf,TM_setup_adj,i_out,n_block,tol);
            f = B*f(1:nq,:); 

            f_hat = fft(f,[],2)*dt_sampling;
            f_hat = invW*f_hat;
            
            % De-normalize frequenci components
            for i=1:nf
                f_hat(:,i)  = f_hat(:,i)*norms(i);
            end

%             subplot(212)
% %                 [f_hat,-L0'\q_hat]
% 
%                 fadj = @(t) f_hat             *exp(2i*pi*fList'*t);
%                 tt = (0:size(fhist,2)-1)*dt;
%                 ftest = fadj(tt );
%                 plot(tt,fhist(1,:),tt,ftest(1,:),':',[i_out+n_block*1]*dt,f(1,:),'o');

            f_out(:,i_m,:,i_iter) = f_hat;

        end
        %% Make next inputs ortogonal to all previous ones
        for i_f = 1:nf
            f  = reshape( f_out(:,:,i_f,i_iter  ),nq,m       ) ;  % current outputs
            F  = reshape( f_in (:,:,i_f,1:i_iter),nq,m*i_iter) ;    %previous input
            f = f - F*(F'*f);  % keep only f that is orthogonal to F 
            [f,~]= qr(f,0);
            f_prev(:,:,i_f) = f;
        end
    end
    
    %% Compute gains based on the Arnold decomposition constructed
    S       = zeros(n_Iter,nf)      ; %gains 
    SS_conv = nan(n_Iter,n_Iter,nf) ; %gains convergence history

    for i_f = 1:nf
        F_in   = reshape( f_in  (:,:,i_f,:),nq,m*n_Iter) ;    % Inputs
        F_out  = reshape( f_out (:,:,i_f,:),nq,m*n_Iter) ;    % Output
        H = F_in'*F_out;
        
        [psi,SS]        = eig(H);
        SS              = sqrt(diag(SS));
        [S(:,i_f),order]= sort(SS,'descend');
        V(:,:,i_f)      = F_in*psi(:,order); 
        
        for j=1:n_Iter
            SS_conv(1:j,j,i_f) = sort(sqrt( eig(H(1:j,1:j))) ,'descend');
        end
    end
    
end

