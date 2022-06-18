function [V,LAMBDA,resHist] = block_eigs(fun,n,neigs,nkirlov,mkirlov,pkirlov,tol)

    np_kirlov=nkirlov*pkirlov; 
    v = zeros(n,np_kirlov);
    w = zeros(n,np_kirlov);
    
    
    %set initial ortogonal random space.
    [fk,~] = qr(randn(n,pkirlov),0);

    i_kirlov = 0 ;  % initialize index for kirlov space 
    resHist=[];
    for i_iter = 1:300
    while i_kirlov < np_kirlov
        n_space  = i_kirlov;
        i_start = n_space+1;
        i_end   = n_space+pkirlov;
        
        i_prev = 1:n_space;
        i_curr = i_start:i_end ;
        i_all  = 1:i_end ;
        
        % add new vector and ortogonalize
        v(:,i_curr)     = fk;
        [v(:,i_all),R]  = qr(v(:,i_all),0);
        w(:,i_prev)     = w(:,i_prev)/R(i_prev,i_prev);

        % Apply linear operator on the new vectors
        for j  = 1:pkirlov
            fprintf('   Filling in w_%.0d\n',i_curr(j))
            w(:,i_curr(j)) = fun(v(:,i_curr(j)));        
        end
        
        
        % project out part of w that in the v subspace                
        n_space  = n_space+pkirlov;

        H  = v(:,1:n_space)' * w(:,1:n_space) ;
        fk = w(:,i_curr)  - v(:,1:n_space)*H(:,end);
        
    
        [psi,lambda] =  eig(H); lambda=diag(lambda);

        [lambda,order] = sort(lambda,'descend');
        psi=psi(:,order);
        V       = v(:,1:n_space)*psi;
        LAMBDA  = lambda;
        
        % Av = w = v*v'*w + (1-v*v')*w
        % Av = w = v*H + f
        % H psi = psi * L
        % A v psi  = v*H*psi + f*psi
        % A V = v*H*psi + f*psi
        % A V = v*psi*L + f*psi
        % A V = V*L + f*psi
        % A V_i = V_i*L_i + (f*psi)_i
        
        i_kirlov = i_kirlov + pkirlov;

        
        f = w-v*v'*w;
        ff = f(:,1:i_kirlov)*psi;
        res = sqrt(sum(abs(ff(:,1:i_kirlov)).^2,1));
        resHist(1:i_kirlov,end+1)=sort(abs(res));
%         norm(fk * psi(end,:))
%         res = abs(norm(fk) * psi(end,:));
%         resHist(1:i_kirlov+1,end+1)=sort(abs(res));
        p=(abs(res(1:i_kirlov))>tol);
        n_conv = sum(~p);
        fprintf('Iter %d , %d convergend modes. Min res %e , min non-converged res %e, |fk|=%e \n',i_iter,n_conv,min(res),min(res(p)),norm(fk(:,end)));
        if (n_conv >= neigs); return; end


    end
    
    %Construct sketch of A
    
    
    
    % Implicit shift and invert
    lambdatar = sort(lambda); 
    
    lambdatar=lambdatar(1:mkirlov*pkirlov);
    
    ilamb = 0;
    Q = eye(size(H));
    for i=1:mkirlov*pkirlov
        l = lambdatar(i);
        H0 = H;
        [Q,R] = qr(H(1:i_kirlov,1:i_kirlov)-eye(i_kirlov)*l);
        
        v(:,1:i_kirlov) = v(:,1:i_kirlov)*Q;
        w(:,1:i_kirlov) = w(:,1:i_kirlov)*Q;
    
        v(:,i_kirlov) = 0 ;
        w(:,i_kirlov) = 0 ;
        i_kirlov = i_kirlov-1;
        
        H  = v(:,1:i_kirlov)'*w(:,1:i_kirlov);
        fk = w(:,i_kirlov-pkirlov+1:i_kirlov);
        
    end 
    end
    