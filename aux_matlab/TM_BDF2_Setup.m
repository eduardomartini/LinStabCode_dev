function [TM_setup,TM_setup_adj,Al,Ar] = TM_BDF2_Setup(H,L,dt,maxIter,tol,toliLU,verbose)
    
    if (~exist('maxIter') && isempty(maxIter)); maxIter   =100   ; end
    if (~exist('tol    ') && isempty(tol    )); tol       =1e-6  ; end
    if (~exist('toliLU ') && isempty(toliLU )); toliLU    =1e-4  ; end
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
    
    % Setup preconditioner for Ar, based on a pivoted iLU decomposition 
    Pi = speye(n,1)*0;
    Pj = speye(n,1)*0;

    pj = symrcm(Al_11 ); % Other tests orderings : symrcm ,  amd,symamd,dissect,colperm,colamd
    pi = symrcm(Al_11'); % 
    
    for i=1:n
        Pi(i,pj(i)) = 1; 
        Pj(pi(i),i) = 1; 
    end

    [iLL,iUU]    = ilu(Al_11(pi,pj),struct('type','ilutp','droptol',toliLU,'milu','off'));

    prec     = @(x) Pi'*(iUU \(iLL \(Pj'*(x))));
    prec_H   = @(x) Pj *(iLL'\(iUU'\(Pi *(x)))); 


    % Function to solve Al*x = y for x
    if verbose==2
        invAl_fun    = @(x) [      cgs(Al_11 ,x(1:n),tol,maxIter,prec  );x(n+1:2*n)];
        invAl_H_fun  = @(x) [      cgs(Al_11',x(1:n),tol,maxIter,prec_H);x(n+1:2*n)];
    else
        invAl_fun    = @(x) [quiet_cgs(Al_11 ,x(1:n),tol,maxIter,prec  );x(n+1:2*n)];
        invAl_H_fun  = @(x) [quiet_cgs(Al_11',x(1:n),tol,maxIter,prec_H);x(n+1:2*n)];
    end
    Ar_fun       = @(x) Ar *x(:) ;
    Ar_H_fun     = @(x) Ar'*x(:) ;
    
    TM_setup.verbose       = verbose; 
    TM_setup.n             = size(L,1);
    TM_setup.n_multi       = 2;
     
    TM_setup.invAl_fun     = invAl_fun;  
    TM_setup.Ar_fun        = Ar_fun;  
    
    TM_setup.dt            = dt;
    TM_setup.setup_ic      = @(x) x(:);
    TM_setup.setup_f       = @(f,t) [(dt*2/3)*f( t+dt );zeros(n,1)];

    TM_setup_adj           = TM_setup;
    TM_setup_adj.invAl_fun = invAl_H_fun;  
    TM_setup_adj.Ar_fun    = Ar_H_fun;  

    disp(   ['Done in ' datestr(toc/24/3600, 'HH:MM:SS')  ' \n'] )

