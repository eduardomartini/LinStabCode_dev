function [q_hat,q0] = CRB_TM_gmres(TM_setup,fList,f_hat,B,C,i_out,n_block,tol,flagAdj)
    TM_setup.verbose = false; % forces quiet run
    
    q0 = zeros(TM_setup.n*TM_setup.n_multi,1); 
    % Compute F
    [q_hat0,qf] = CRB_TM(TM_setup,fList,q0,f_hat,B,C,i_out,n_block,1,-inf,flagAdj);
    % Solve Ax = F
    % qf = exp(-At)q0 + F
    % (1-exp(-At))q0 -F = (q0-qf) -F = fk 
    % GMRES
    F = -(q0-qf) ; % get RHS of the linear problem.
    fk = F ; %initialize residual
    err = inf; %
    
    Q=[]; %initialize snapshot matrix
    i=0;
    W=[];
    while err>tol
        [Q,R]=qr([Q,fk],0); % assures that components are ortonormal;
        W=W/R(1:end-1,1:end-1); 
        
        q0 = Q(:,end);
        
        i = i + 1;
%         [q_hats{i},qf] = CRB_TM(TM_setup,fList,q0,f_hat,B,C,i_out,n_block,1,-inf,flagAdj);
%         qf=qf-F;
        [q_hats{i},qf] = CRB_TM(TM_setup,fList,q0,f_hat*0,B,C,i_out,n_block,1,-inf,flagAdj);
        qf=qf;
        fk = (q0-qf);
        W(:,end+1)=fk; 
        
        psi = pinv(W) * F;
        err = (norm(W*psi-F)/norm(F));
    end
    disp(['GMRES error at iteration ' num2str(i) ' : ' num2str(err,'%3.1e' )]);
    % Reconstruct solution
    q_hat = q_hat0;
    for i=1:length(psi)
        q_hat = q_hat + q_hats{i}*psi(i);
    end
    
    %% Get fourrier components (partial solutions are not periodic, so we cant trust FFT there..)
    % reconstruct initial condition
%     q0    = Q*psi;  
%     [q_hat2,qf] = CRB_TM(TM_setup,fList,q0,f_hat,B,C,i_out,n_block,1,-inf,flagAdj);
%     norm(q0-qf) % check periodicity assumption
%     for i=1:length(fList)
%         norm(q_hat(:,i)-q_hat2(:,i))
%     end

%     [q_hat2,qf2] = CRB_TM(TM_setup,fList,q0*0,f_hat,B,C,i_out,n_block,inf,tol,flagAdj);
%     [q_hat3,qf3] = CRB_TM(TM_setup,fList,qf2,f_hat,B,C,i_out,n_block,inf,tol,flagAdj);
% 
%     [q_hat4,qf4] = CRB_TM(TM_setup,fList,qf2,f_hat,B,C,i_out,n_block,1,-inf,flagAdj);
%     [q_hat5,qf5] = CRB_TM(TM_setup,fList,qf4,f_hat,B,C,i_out,n_block,1,-inf,flagAdj);

%     q_hat = q_hat2;
        
    end