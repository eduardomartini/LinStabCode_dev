function [TM_setup,TM_setup_adj,Al,Ar] = TM_EulerImplicit_Setup(H,L,dt,verbose,opts,filters)
    
    if ~exist('verbose','var') ; verbose   =false ; end
    if ~exist('opts','var')    ; opts=struct('type','ilu');end
%     if ~exist('maxIter','var') ; maxIter   =100   ; end
%     if ~exist('tol    ','var') ; tol       =1e-6  ; end
%     if ~exist('toliLU ','var') ; toliLU    =1e-4  ; end
    fprintf('Setting up backwards Euler integration scheme...');
    tic;
    n = size(L,1);
    % Setup Al and Ar matrices
    Al = ( H - dt * L );
      
    Ar = ( H  );
 
    [invAl_fun,invAl_H_fun] = GetInverseFunction(Al,opts);

    Ar_fun       = @(x) Ar *x;
    Ar_H_fun     = @(x) Ar'*x;
    
    TM_setup.verbose        = verbose; 
    TM_setup.n              = size(L,1);
    TM_setup.n_multi        = 1;
    
    if exist('filters','var') 
        FILTER    = @(x) reshape( filters.filter   ( reshape(x,[],5)),[],1) ; 
        FILTER_ct = @(x) reshape( filters.filter_ct( reshape(x,[],5)),[],1) ; 
    else
        FILTER    = @(x) x ; 
        FILTER_ct = @(x) x ; 
    end


    TM_setup.invAl_fun      = @(x) FILTER(invAl_fun(x)) ;
    TM_setup.Ar_fun         = Ar_fun;  
    
    TM_setup.dt             = dt;
    
    TM_setup.setup_ic       = @( x ) x;
    TM_setup.setup_f        = @(f,t) dt*f(t+dt);
    

    TM_setup_adj            = TM_setup;
    TM_setup_adj.invAl_fun  = @(x) FILTER_ct(invAl_H_fun(x));  
    
    TM_setup_adj.Ar_fun     = Ar_H_fun;

    disp(   ['Done in' datestr(toc/24/3600, 'HH:MM:SS')  ' \n'] )

