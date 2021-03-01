function [err,d_inf] = Solution_Statistics(A_t,y,u,p,y_hat,tau,hx,hy,gamma,ny,nx,nt,obj_norm,M1,scale_const,w_y,w_u);
% ============================================================================================================== %
% This function takes as an input the solution and the data of the test problem,
% and returns the statistics of the solution:
%                  1) err; the discrepancy e(y-\bar{y}), measured in the norm indicated by (obj_norm),
%                  2) d_inf; the dual infeasibility.
% ______________________________________________________________________________________________________________ %
    if (nargin < 16 || isempty(w_y))
        w_y = zeros(ny*nx*nt,1);
    end
    if (nargin < 17 || isempty(w_u))
        w_u = zeros(ny*nx*nt,1); 
    end
    lvls = [ny;nx;nt];
    A_t = reshape(A_t,ny*nx*(2*nt-1),1);
    A_tcmat = Multilevel_Circulant_Extrapolation(A_t, lvls, 3);
    err = norm((y(1:nx*ny*(nt-1),1) - y_hat(1:nx*ny*(nt-1),1)))^2;
    if (obj_norm ~= 1)
        
        %% Approximation of the state discrepancy (L^2 norm) %%
        %err = err + 0.5 * norm((y(nx*ny*(nt-1)+1:end,1) - (2.*y_hat(nx*ny*(nt-1)+1:end,1))))^2;
                err = err + 0.5 * norm((y(nx*ny*(nt-1)+1:end,1) - (y_hat(nx*ny*(nt-1)+1:end,1))))^2;
        err = sqrt(tau*hx*hy*err);
        %err_scale = sqrt(norm(y_hat(1:nx*ny*(nt-1),1))^2 + norm(2.*y_hat(nx*ny*(nt-1)+1:end,1))^2);
        %err = err/err_scale;
        
        %% Computation of dual infeasibility (L^2 norm) 
       % if (obj_norm ~= 0)
            tmp_vec =  (y-y_hat);  
            tmp_vec(nx*ny*(nt-1)+1:end,1) = (1/2).*(y(nx*ny*(nt-1)+1:end,1) - 2.*y_hat(nx*ny*(nt-1)+1:end,1));
            d_inf = norm(Lvl3_Toeplitz_Operator(A_tcmat,p,ny,nx,nt) + tmp_vec + w_y,Inf);
            tmp_vec = gamma.*u;
            tmp_vec(nx*ny*(nt-1)+1:end,1) = (1/2).*tmp_vec(nx*ny*(nt-1)+1:end,1);
            d_inf = max(d_inf,norm(tmp_vec + scale_const.*p + scale_const.*w_u,Inf));
     %   end
    else
        %% Approximation of the state discrepancy (H^1 norm)
        err = err + 0.5 * norm(y(nx*ny*(nt-1)+1:end,1) - y_hat(nx*ny*(nt-1)+1:end,1));
        err = sqrt(tau^2*err);
        M1 = reshape(M1,ny*nx*(2*nt-1),1);
        M1_cmat = Multilevel_Circulant_Extrapolation(M1, lvls, 3);
        tmp_2 = y-y_hat;
        tmp_1 = Lvl3_Toeplitz_Operator(M1_cmat,tmp_2,ny,nx,nt);
        tmp_1(nx*ny*(nt-1)+1:end,1) = 0.5 .* tmp_1(nx*ny*(nt-1)+1:end,1);
        tmp_err = tau^2*(tmp_2'*tmp_1);
        err = err + tmp_err;
        
        
        %% Computation of dual infeasibility (H^1 norm)
        
        d_inf = norm(Lvl3_Toeplitz_Operator(A_tcmat,p,ny,nx,nt) + tmp_1,Inf);
        tmp_vec = gamma.*u;
        tmp_vec(nx*ny*(nt-1)+1:end,1) = (1/2).*tmp_vec(nx*ny*(nt-1)+1:end,1);
        d_inf = max(d_inf,norm(tmp_vec + scale_const.*p,Inf));
        
        
    end
      %     d_inf = 0;

end

