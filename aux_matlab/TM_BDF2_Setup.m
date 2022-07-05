function [TM_setup,TM_setup_adj,Al,Ar] = TM_BDF2_Setup(H,L,dt,verbose,opts,filters)
    
    if ~exist('opts','var')    ; opts=struct('type','ilu');end
    if (~exist('verbose') && isempty(verbose)); verbose   = 1    ; end
    fprintf('Setting up 2-nd order backwards diferentiation scheme (BDF2)...');
    tic;
    n = size(L,1);
    
    %% Differential equation : H dq/dt = A q + f
    %% Integration sceheme   : H (q_n+1-4/3*q_n+1/3*q_n-1)/dt = 2/3*(A*q_n+1 f_n+1) 
    %% Iteration             : [q_n+1 ; q_n] = Ar^-1 * (Al*[q_n;q_n-1] + [f_n+1;0])
    
    % Setup Al and Ar matrices
    I = speye(n);
    Z =     I*0 ;

    %Al    = [ H - dt * L*2/3 , Z ; Z , I ] );
    Al_11 = [ H-(dt*2/3)*L ];
    Ar    = [ (4/3)*H , (-1/3)*H ; I , Z ];
    
    
    [invAl11_fun,invAl11_H_fun] = GetInverseFunction(Al_11,opts);


    if exist('filters','var') 
        FILTER    = @(x) reshape( filters.filter   ( reshape(x,[],5)),[],1) ; 
        FILTER_ct = @(x) reshape( filters.filter_ct( reshape(x,[],5)),[],1) ; 
    else
        FILTER    = @(x) x ; 
        FILTER_ct = @(x) x ; 
    end

    invAl_fun       = @(x) [ FILTER(   invAl11_fun  (x(1:n))) ; x(n+1:2*n)];
    invAl_ct_fun    = @(x) [ FILTER_ct(invAl11_H_fun(x(1:n))) ; x(n+1:2*n)];


    Ar_fun       = @(x) Ar *x(:) ;
    Ar_ct_fun    = @(x) Ar'*x(:) ;
    
    TM_setup.verbose       = verbose; 
    TM_setup.n             = size(L,1);
    TM_setup.n_multi       = 2;
     
    TM_setup.invAl_fun     = invAl_fun;  
    TM_setup.Ar_fun        = Ar_fun;  
    
    TM_setup.dt            = dt;
    TM_setup.setup_ic      = @(x) x(:);
    TM_setup.setup_f       = @(f,t) [(dt*2/3)*f( t+dt );zeros(n,1)];

    TM_setup_adj           = TM_setup;
    TM_setup_adj.invAl_fun = invAl_ct_fun;  
    TM_setup_adj.Ar_fun    = Ar_ct_fun;  

    disp(   ['Done in ' datestr(toc/24/3600, 'HH:MM:SS')  ' \n'] )

