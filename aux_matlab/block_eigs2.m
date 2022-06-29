function [V,LAMBDA_conv,LAMBDA_all,res,resHist] = block_eigs2(applyOP,n,neigs,nkirlov,mkirlov,pkirlov,tol)
%     umf = matlab.internal.decomposition.builtin.UMFPACKWrapper(A, 0.1, 0.001, true, 0);

    np_kirlov=round(nkirlov/pkirlov)*pkirlov; 
    mkirlov=round(mkirlov/pkirlov)*pkirlov; 
    v = zeros(n,np_kirlov);
    w = zeros(n,np_kirlov);
    
    H = zeros(np_kirlov);
    
    %set initial ortogonal random space.
    [fk,~] = qr(randn(n,pkirlov),0);

    i_kirlov = 0 ;  % initialize index for kirlov space 
    resHist=[];
    for i_iter = 1:300
    while i_kirlov < np_kirlov
        i_start = i_kirlov+1;
        i_end   = i_kirlov+pkirlov;
        
        i_prev = 1:i_kirlov;
        i_curr = i_start:i_end ;
        i_all  = 1:i_end ;
        
        % add new vector and ortogonalize
%         v(:,i_curr)     = qr(fk-v(:,i_all)*(pinv(v(:,i_all))*fk),0);

%         v(:,i_curr)     = fk;
%         for i=i_curr
%             [projv, dw] = matlab.internal.math.projectInto(v, v(:,i), i);
%             v(:,i) = v(:,i)-projv
%         end

%         v(:,i_curr)     = fk;
%         for i=i_curr
%            [proj2 , ~] = matlab.internal.math.projectInto(v, v(:,i), i-1);
%             v(:,i) = v(:,i) - proj2; 
%             v(:,i)=v(:,i)/norm(v(:,i));
%         end

        [v(:,i_curr),R]     = qr(fk,0);

        [v(:,i_all),R]     = qr(v(:,i_all),0);
        w(:,i_prev)     = w(:,i_prev)/R(i_prev,i_prev);
        
        input  = v(:,i_curr);
        output = zeros(size(input));
        if pkirlov==1
            for j  = 1:pkirlov
               output(:,j) = applyOP(input(:,j));        
            end
            w(:,i_curr) = output;
        else    
            for j  = 1:pkirlov
               output(:,j) = applyOP(input(:,j));        
            end
            w(:,i_curr) = output;
        end


        % project out part of w that in the v subspace                
%         for i=i_curr
%             for j=1:i
%                 H(j,i)  = v(:,j)' * w(:,i) ;
%             end
%             for j=1:min(pkirlov,i-1) 
%                 H(i,i-j)  = v(:,i)' * w(:,i-j) ;
%             end
%         end
% 
%         for j=1:pkirlov
%             [proj , ~] = matlab.internal.math.projectInto(v, w(:,i_curr(j)), i_kirlov+pkirlov);
%             fk(:,j) = w(:,i_curr(j)) - proj;
%         end

        i_kirlov = i_kirlov + pkirlov;

        H(1:i_kirlov,1:i_kirlov)  = v(:,1:i_kirlov)' * w(:,1:i_kirlov) ;
        fk = w(:,i_curr)  - v(:,1:i_kirlov)*H(1:i_kirlov,i_kirlov);

        check_conv = mod(i_kirlov,1)==0;
        if check_conv
            [psi,lambda] =  eig(H(1:i_kirlov,1:i_kirlov)); lambda=diag(lambda);
            [lambda,order] = sort(lambda,'descend');
            psi=psi(:,order);

            % Av = w = v*v'*w + (1-v*v')*w
            % Av = w = v*H + f
            % H psi = psi * L
            % A v psi  = v*H*psi + f*psi
            % A V = v*H*psi + f*psi
            % A V = v*psi*L + f*psi
            % A V = V*L + f*psi
            % A V_i = V_i*L_i + (f*psi)_i




        % Compute full residualds
    %         f = w-v*v'*w;
    %         ff = f(:,1:i_kirlov)*psi;
        % Compute only relevant residualds
            f =  fk;
            ff = f*psi(end-pkirlov+1:end,:);

            res = zeros(i_kirlov,1);
            for i=1:i_kirlov
                res(i) = norm(ff(:,i));
            end
            

    %         res = sqrt(sum(abs(ff(:,1:i_kirlov)).^2,1));
            resHist(1:i_kirlov,end+1)=sort(abs(res));
    %         norm(fk * psi(end,:))
    %         res = abs(norm(fk) * psi(end,:));
    %         resHist(1:i_kirlov+1,end+1)=sort(abs(res));
            p=(abs(res(1:i_kirlov))>tol);
            n_conv = sum(~p);
            fprintf('Iter %d, %d , %d convergend modes. Min res %e , min non-converged res %e, |fk|=%e \n',i_iter,i_curr(end),n_conv,min(res),min(res(p)),norm(fk(:,end)));
            if (n_conv >= neigs) 
                V       = v(:,1:i_kirlov)*psi;
                LAMBDA  = lambda;

                LAMBDA_conv = LAMBDA(1:neigs);
                LAMBDA_all  = LAMBDA;
                return; 
            end
        end

    end
    
    % Implicit shift and invert
    lambdatar = sort(lambda); 
    
    lambdatar=lambdatar(1:mkirlov);
    
    ilamb = 0;
    Q = eye(size(H));
    for i=1:mkirlov
        l = lambdatar(i);
        H0 = H;
        [Q,R] = qr(H(1:i_kirlov,1:i_kirlov)-eye(i_kirlov)*l);
        
        v(:,1:i_kirlov) = v(:,1:i_kirlov)*Q;
        w(:,1:i_kirlov) = w(:,1:i_kirlov)*Q;
    
        v(:,i_kirlov) = 0 ;
        w(:,i_kirlov) = 0 ;
        i_kirlov = i_kirlov-1;
        
        H(1:i_kirlov,1:i_kirlov)  = v(:,1:i_kirlov)'*w(:,1:i_kirlov);        
    end 
        
    fk = w(:,i_kirlov-pkirlov+1:i_kirlov);
    for j=1:pkirlov
        [proj , ~] = matlab.internal.math.projectInto(v, w(:,i_kirlov-pkirlov+1), i_kirlov);
        fk(:,j) =fk(:,j) - proj;
    end
%     fk = w(:,i_kirlov-pkirlov+1:i_kirlov);
%     [psi,lambda] =  eig(H(1:n_space,1:n_space)); lambda=diag(lambda);

    end
    