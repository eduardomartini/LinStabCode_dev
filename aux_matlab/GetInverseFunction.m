function [invA_fun,invA_T_fun] = GetInverseFunction(A,opts)%,tol,toliLU,maxIter,verbose)
    if ~exist('opts','var' ); opts=struct(); end
    if ~isfield(opts,'type'); opts.type = 'lu'; end
        
    if strcmp(opts.type,'builtin')
        disp('Computing inverse function using matlab internal function')
        umf         = matlab.internal.decomposition.builtin.UMFPACKWrapper(A , 0.1, 0.001, true, 0);
        umf_T       = matlab.internal.decomposition.builtin.UMFPACKWrapper(A', 0.1, 0.001, true, 0);
        invA_fun    = @(v) solve(umf, v, false);
        invA_T_fun  = @(v) solve(umf_T, v, false);
    elseif strcmp(opts.type,'lu')
        disp('Computing inverse function LU decomposition')
        [LL,UU,pp,qq,rr]    = lu(A);
        invA_fun   = @(v) qq *(UU \(LL \(pp *(rr \v))));
        invA_T_fun = @(v) rr'\(pp'*(LL'\(UU'\(qq'*v))));
    elseif strcmp(opts.type,'ilu')
        disp('Computing inverse function iLU preconditioner and an iterative method;')
        if ~isfield(opts,'maxIter') ; opts.maxIter   =100    ; end
        if ~isfield(opts,'tol'    ) ; opts.tol       =1e-6   ; end
        if ~isfield(opts,'toliLU ') ; opts.toliLU    =1e-4   ; end
        if ~isfield(opts,'verbose') ; opts.verbose   =false  ; end
        if ~isfield(opts,'solver' ) ; opts.solver    ='gmres'; end


        % Setup preconditioner for Ar
        n=size(A,1);
        %Create pivoting matrices/indexes
        Pi = speye(n,1)*0;
        Pj = speye(n,1)*0;

        pj = symrcm(A ); % symrcm ,  amd,symamd,dissect,colperm,colamd
        pi = symrcm(A'); % symrcm ,  amd,symamd
        for i=1:n
            Pi(i,pj(i)) = 1; 
            Pj(pi(i),i) = 1; 
        end

        % Get iLU decomposition to be used as preconditioner
        [iLL,iUU]    = ilu(A(pi,pj),struct('type','ilutp','droptol',opts.toliLU,'milu','off'));

        prec     = @(x) Pi'*(iUU \(iLL \(Pj'*(x))));
        prec_H   = @(x) Pj *(iLL'\(iUU'\(Pi *(x)))); 

        % Function to solve Al*x = y for x
        if opts.verbose
            if strcmp('cgs',opts.solver)
                invA_fun    = @(x)       cgs(A  ,x,opts.tol,opts.maxIter,prec);
                invA_T_fun  = @(x)       cgs(A' ,x,opts.tol,opts.maxIter,prec_H);
            elseif strcmp('gmres',opts.solver)
                invA_fun    = @(x)       gmres(A  ,x,[],opts.tol,opts.maxIter,prec);
                invA_T_fun  = @(x)       gmres(A' ,x,[],opts.tol,opts.maxIter,prec_H);
            end
        else
            if strcmp('cgs',opts.solver)
                invA_fun    = @(x) quiet_cgs(A ,x,opts.tol,opts.maxIter,prec);
                invA_T_fun  = @(x) quiet_cgs(A',x,opts.tol,opts.maxIter,prec_H);
            elseif strcmp('gmres',opts.solver)
                invA_fun    = @(x) quiet_gmres(A  ,x,[],opts.tol,opts.maxIter,prec);
                invA_T_fun  = @(x) quiet_gmres(A' ,x,[],opts.tol,opts.maxIter,prec_H);
            end
        end
    end