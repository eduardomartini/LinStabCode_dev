function [q_out,varargout] = TM(q0,f,nPeriods,TM_setup,i_block_output,n_block,tol,filter_fun)
    flag_saveHist =  nargout==2 ; 

    if ~exist('filter_fun'); filter_fun = @(x) x;end
    Ar      = TM_setup.Ar_fun;
    invAl   = TM_setup.invAl_fun;
    dt      = TM_setup.dt;
    verbose = TM_setup.verbose;
    n        =TM_setup.n;

    err = 1e99;
    i_iter=0;
    err = inf;
    
    if verbose==2; str_end='\n ';
    else         ; str_end='';
    end
    
    if verbose>0 
        slen=fprintf('Current time(iter) %3.1e(%.0f), log(err) %2.1e %s ',i_iter*dt,i_iter,err,str_end) ;
    end
    
    if flag_saveHist
        if isinf(nPeriods) ; qhist=nan(n,1);
        else               ; qhist=nan(n,nPeriods*n_block);
        end
    end
    
    n_outputs = length(i_block_output);
    q_out = zeros(length(q0),n_outputs);
    
    q       = TM_setup.setup_ic(q0);
    f_fun   = TM_setup.setup_f;
    
    q_check    = q ; % vector containing previous state for error estimation
    q_out(:,1) = q ; % initialize output vector
    while err>tol && i_iter/n_block <nPeriods
        time = i_iter*dt;
        i_iter=i_iter+1;
        
        qold = q;

        q = invAl ( Ar(q)+f_fun(f,time));
        q = filter_fun(q);
        
        if mod(i_iter,n_block)==0
            err = norm(q_check-q)/max(norm(q),norm(q_check));
            q_check=q;
        end
        
        i_block = mod(i_iter,n_block); % time step position in a block
        iout = find( i_block==i_block_output);
        if ~isempty(iout)
            q_out(:,iout) = q;
        end

            
        if verbose==1; fprintf(repmat('\b',1,slen)) ;end
        if verbose>0 
            slen=fprintf('Current time(iter) %3.1e(%.0f), log(err) %2.1e %s ',i_iter*dt,i_iter,err,str_end) ;
        end
        
        if flag_saveHist
            qhist(:,i_iter+1) = q(1:n);
        end
    end
    
    if flag_saveHist
        varargout{1}=qhist;
    end
    
    if verbose>0
        disp(' Done!')
    end