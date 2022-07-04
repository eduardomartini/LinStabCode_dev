function [q_hat,qf] = CRB_TM(TM_setup,fList,q0,f_prev,B,C,i_out,n_block,nPeriods,tol,flagAdj)
        nf=length(fList);        
        if nf==1
            dt_sampling = 1;
        else
            i_sampling = i_out(2);
            dt_sampling = TM_setup.dt*i_sampling;
        end
        
        if flagAdj
            f = @(t) (B*f_prev)*exp(-2i*pi*fList.'*t);
        else
            f = @(t) (C*f_prev)*exp( 2i*pi*fList.'*t);
        end
        
        q = TM(q0,f ,nPeriods,TM_setup,i_out,n_block,tol);
        qf = q(:,1);
        %%
        if flagAdj==false
            q     = C*q(1:TM_setup.n,:);
            q_hat = fft(q,[],2)/nf;            
        else
            q     = B*q(1:TM_setup.n,:);
            q_hat = ifft(q,[],2);%/nf;            
        end
        
        
%         q_hat = q_hat;
            
end