function [TM_setup,TM_setup_adj,Al,Ar] = TM_EulerImplicit_Setup(H,L,dt,maxIter,tol,toliLU,verbose)
    
    if ~exist('maxIter') ; maxIter   =100   ; end
    if ~exist('tol    ') ; tol       =1e-6  ; end
    if ~exist('toliLU ') ; toliLU    =1e-4  ; end
    if ~exist('verbose') ; verbose   =false ; end
    fprintf('Setting up backwards Euler integration scheme...');
    tic;
    n = size(L,1);
    % Setup Al and Ar matrices
    Al = ( H - dt * L );
    Ar = ( H  );
    
    

    % Setup preconditioner for Ar, based on a pivoted iLU decomposition 
    Pi = speye(n,1)*0;
    Pj = speye(n,1)*0;

    pj = symrcm(Al ); % symrcm ,  amd,symamd,dissect,colperm,colamd
    pi = symrcm(Al'); % symrcm ,  amd,symamd
    for i=1:n
        Pi(i,pj(i)) = 1; 
        Pj(pi(i),i) = 1; 
    end

    [iLL,iUU]    = ilu(Al(pi,pj),struct('type','ilutp','droptol',toliLU,'milu','off'));

    prec     = @(x) Pi'*(iUU \(iLL \(Pj'*(x))));
    prec_H   = @(x) Pj *(iLL'\(iUU'\(Pi *(x)))); 

    % Function to solve Al*x = y for x
    if verbose == 2
        invAl_fun    = @(x)       cgs(Al  ,x,tol,maxIter,prec);
        invAl_H_fun  = @(x)       cgs(Al' ,x,tol,maxIter,prec_H);
    else
        invAl_fun    = @(x) quiet_cgs(Al ,x,tol,maxIter,prec);
        invAl_H_fun  = @(x) quiet_cgs(Al' ,x,tol,maxIter,prec_H);
    end
    
    Ar_fun       = @(x) Ar *x;
    Ar_H_fun     = @(x) Ar'*x;
    
    TM_setup.verbose   = verbose; 
    TM_setup.n         = size(L,1);
    TM_setup.n_multi   = 1;
    
    TM_setup.invAl_fun   = invAl_fun;  
    TM_setup.Ar_fun      = Ar_fun;  
    
    TM_setup.dt        = dt;
    
    TM_setup.setup_ic  = @( x ) x;
    TM_setup.setup_f   = @(f,t) dt*f(t+dt);
    

    TM_setup_adj = TM_setup;
    TM_setup_adj.invAl_fun = invAl_H_fun;  
    TM_setup_adj.Ar_fun    = Ar_H_fun;  

    disp(   ['Done in' datestr(toc/24/3600, 'HH:MM:SS')  ' \n'] )

