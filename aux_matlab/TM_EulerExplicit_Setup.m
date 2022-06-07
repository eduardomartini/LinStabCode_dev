function [TM_setup,Al_M,Ar_M] = TM_EulerImplicit_Setup(H,L,dt,maxIter,tol,toliLU,verbose)
    
    if ~exist('maxIter') ; maxIter   =100   ; end
    if ~exist('tol    ') ; tol       =1e-6  ; end
    if ~exist('toliLU ') ; toliLU    =1e-4  ; end
    if ~exist('verbose') ; verbose   =false ; end
    fprintf('Setting up backwards Euler integration scheme...');
    tic;
    n = size(L,1);
    % Setup Al and Ar matrices
    
    Ar_fun       = @(x) (x+L *x*dt);
    Ar_H_fun     = @(x) (x+L'*x*dt);
    invAl_fun    = @(x) x; 
    
    TM_setup.verbose   = verbose;     
    TM_setup.n         = size(L,1);
    TM_setup.n_multi   = 2;

    TM_setup.invAl_fun = invAl_fun;  
    TM_setup.Ar_fun    = Ar_fun;  
    
    TM_setup.dt        = dt;
    TM_setup.setup_ic  = @( x ) x;
    TM_setup.setup_f   = @(f,t) dt*f(t);
    
    
    TM_setup_adj           = TM_setup;
    TM_setup_adj.invAl_fun = invAl_fun;  
    TM_setup_adj.Ar_fun    = Ar_H_fun;  

    
    disp(   ['Done in' datestr(toc/24/3600, 'HH:MM:SS')  ' \n'] )

