function [y,u,p,w_y,w_u,opt,k,avg_iter] = Box_Constrained_ADMM(A,A_t,g,y_hat,gamma,ny,nx,nt,P,mode,lb,ub,delta,scale_const,beta)
% ==================================================================================== %
% [x,y] = P_ADMM(A,Q,b,lb,ub): This is a Proximal-Alternating Direction Method of 
%                              multipliers, suitable for convex quadratic programming 
%                              problems with linear constraints.
% ------------------------------------------------------------------------------------ %
% Input: A, the constraint matrix,
%        Q, the Hessian,
%        c, the linear part of the objective,
%        b, the right hand side of the equality constraints,
%        lb, ub, vectors containing the inequality constraint bounds, i.e.
%                             x in [lb,ub] or x(i) in [lb(i),ub(i)].
%        tol (default = 10^(-4)), is the tolerance to which the problem is solved.
% Output: x, the optimal solution of the problem,
%         y, the optimal Lagrange multiplier vector.
% ____________________________________________________________________________________ %
    lvls = [ny;nx;nt];
    A = reshape(A,ny*nx*(2*nt-1),1);
    A_t = reshape(A_t,ny*nx*(2*nt-1),1);
    Acmat = Multilevel_Circulant_Extrapolation(A, lvls, 3);
    A_tcmat = Multilevel_Circulant_Extrapolation(A_t, lvls, 3);
    N = nx*ny*nt;  
    % =================================================================================== %
    % Starting point.
    % ----------------------------------------------------------------------------------- %
    y = zeros(N,1);
    u = zeros(N,1);
    if (mode == 2 || mode == 3)
        z_y = zeros(N,1);
        w_y = zeros(N,1);
    else
        w_y = [];
    end
    if (mode == 3 || mode == 4)
        z_u = zeros(N,1);
        w_u = zeros(N,1);
    else
        w_u = [];
    end
    p = zeros(N,1);
    % ___________________________________________________________________________________ %
    header();
    out_iter = 1000;
    in_iter = 300; 
    k = 1;
    out_tol = 1e-4;
    in_tol = 1e-1;
    avg_iter = 0;
    while (k <= out_iter)
   
    % =================================================================================== %
    % 1st ADMM sub-problem: (Solved Inexactly)
    % (x^{k+1},y^{k+1}) = argmin{c^T + (1/2)x^TQx - y^T(Ax-b)-z_k^T(x(cvar)-x_ccp_k)
    %                            + (1/(2*delta))(||x(cvar)-x_ccp_k||^2 + ||Ax-b||^2)}
    % We apply one iteration of the ALM algorithm for this step, i.e. we solve:
    % 
    %       [-(Q+(1/delta)I)   A^T  ] [x]    [c + z_k - (1/delta) x_ccp_k]
    %       [                       ]     = 
    %       [       A        delta I] [y]    [      b + delta y_k        ]
    %
    % ----------------------------------------------------------------------------------- %  
        if (mode == 2)
            tmp_mat = ((scale_const)^2/(beta*gamma)).*ones(N,1); 
            tmp_mat((nt-1)*nx*ny+1:end,1) = ((scale_const)^2/(beta*(gamma/2))).*ones(nx*ny,1); 
            tmp_mat_2 = tmp_mat./scale_const;
            tmp_mat = 1./(tmp_mat + delta/beta);
            rhs = scale_const.*(beta -1)*p;
            rhs = tmp_mat_2.*rhs;
            rhs = rhs + scale_const.*g - delta/beta .* p;  
            rhs = tmp_mat.*rhs;
            rhs = rhs + (1-beta)*p;
            rhs = Lvl3_Toeplitz_Operator(A_tcmat,rhs,ny,nx,nt) + beta*(y_hat - w_y + (1/delta).*z_y);
        elseif (mode == 3 || mode == 4)
            tmp_mat = ((scale_const)^2/(beta*(gamma+(scale_const)^2/delta))).*ones(N,1);
            tmp_mat((nt-1)*nx*ny+1:end,1) = ((scale_const)^2/(beta*(gamma/2+(scale_const)^2/delta))).*ones(nx*ny,1);
            tmp_mat_2 = tmp_mat./scale_const;
            tmp_mat = 1./(tmp_mat + delta/beta); 
            rhs = beta*(scale_const.*w_u - ((scale_const)^2/delta).*z_u) + scale_const.*(beta -1)*p; 
            rhs = tmp_mat_2.*rhs;
            rhs = rhs + scale_const.*g - delta/beta .* p; 
            rhs = tmp_mat.*rhs;
            rhs = rhs + (1-beta)*p;
            if (mode == 3)
                rhs = Lvl3_Toeplitz_Operator(A_tcmat,rhs,ny,nx,nt) + beta*(y_hat - w_y + (1/delta).*z_y);
            else
                rhs = Lvl3_Toeplitz_Operator(A_tcmat,rhs,ny,nx,nt) + beta*y_hat;
            end
        end
        [y,~,~,iter] = pcg(@(w) NE_multiplier(w), rhs, in_tol, in_iter,@(w) Apply_DB_Preconditioner(w));
        avg_iter = avg_iter + iter;
        
        if (mode == 2)
            p_prev = p;
            tmp_vec = Lvl3_Toeplitz_Operator(Acmat,y,ny,nx,nt) - scale_const.*g + (delta/beta).*p - tmp_mat_2.*(scale_const.*(beta-1)*p);
            p = tmp_mat.*tmp_vec;
            u = (tmp_mat_2./scale_const).*(-scale_const.*p + scale_const.*(1-beta)*p_prev);
        elseif (mode == 3 || mode == 4)
            p_prev = p;
            tmp_vec = Lvl3_Toeplitz_Operator(Acmat,y,ny,nx,nt) - scale_const.*g + (delta/beta).*p ...
                      - tmp_mat_2.*(beta*(scale_const.*w_u - ((scale_const)^2/delta).*z_u) + scale_const.*(beta-1)*p); 
            p = tmp_mat.*tmp_vec;
            u = (tmp_mat_2./scale_const).*(-scale_const.*p-beta*(scale_const.*w_u-((scale_const)^2/delta).*z_u) + scale_const.*(1-beta)*p_prev); 
        end
    % ___________________________________________________________________________________ %

    % =================================================================================== %
    % 2nd ADMM sub-problem: (Solved Exactly)
    % (x_ccp^{k+1}) = argmin{ c^Tx_{k+1} + (1/2)x_{k+1}^TQx_{k+1} - y_{k+1}^T(Ax_{k+1}-b)
    %                        -z_k^T(x(cvar)_{k+1}-x_ccp) +
    %                        (1/(2*delta))(||x(cvar)_{k+1}-x_ccp||^2 + ||Ax_{k+1}-b||^2)}
    % ----------------------------------------------------------------------------------- %
        if (mode ~= 4)
            z_y = y + (delta) .* w_y; 
        end
        if (mode == 3 || mode == 4)
            z_u = u + (delta/scale_const).* w_u;
        end
        for i = 1:N  % Stay feasible with respect to box constraints.
            if (mode ~= 4)
                z_y(i) = max((z_y(i)),lb(i,1));
                z_y(i) = min((z_y(i)),ub(i,1));
            end
            if (mode == 3 || mode == 4)
                if (mode == 3)
                    j = N+i;
                else
                    j = i;
                end
                z_u(i) = max((z_u(i)),lb(j,1));
                z_u(i) = min((z_u(i)),ub(j,1));
            end
        end
        
    % ___________________________________________________________________________________ %

    % =================================================================================== %
    % Dual Update:
    % ----------------------------------------------------------------------------------- %
        if (mode ~= 4)
            w_y = w_y + (beta/delta).*(y-z_y); % CHANGE HAPPENED HERE.
        end
        if (mode == 3 || mode == 4)
            w_u = w_u + ((beta*scale_const)/delta).*(u-z_u);
        end
    % ___________________________________________________________________________________ %
    
    % =================================================================================== %
    % Calculate residuals and re-iterate.
    % ----------------------------------------------------------------------------------- %
        res_p_eq = norm(Lvl3_Toeplitz_Operator(Acmat,y,ny,nx,nt)+scale_const.*(u-g),Inf);
        if (mode == 2)
            res_p_ineq = norm((y-z_y),Inf);
        elseif (mode == 3)
            res_p_ineq = (norm(scale_const*(y-z_y),Inf)+norm(scale_const*(u-z_u),Inf));
        elseif (mode == 4)
            res_p_ineq = (norm(scale_const*(u-z_u),Inf));
        end
        
        output(k,res_p_eq,res_p_ineq);
        k = k+1;
        if (res_p_eq <= out_tol && res_p_ineq <= out_tol && res_p_ineq ~= 0)
            opt = 1;
            avg_iter = round(avg_iter/k);
            return;
        elseif (res_p_eq <= out_tol*1e-0 && res_p_ineq <= out_tol)
            opt = 1;
            avg_iter = round(avg_iter/k);
            return;
        end
        if (res_p_ineq > 0)
            in_tol = min(0.05*min(res_p_ineq,res_p_eq),in_tol);
            in_tol = max(in_tol, 1e-6);
        else
            in_tol = min(0.05*res_p_eq,in_tol);
            in_tol = max(in_tol, 1e-6);
        end
    % __________________________________________________________________________________ %
    end
    opt = 0;
    avg_iter = round(avg_iter/k);
        function x = NE_multiplier(w)
            x = Lvl3_Toeplitz_Operator(Acmat,w,ny,nx,nt); 
            x = tmp_mat.*x;
            x = Lvl3_Toeplitz_Operator(A_tcmat,x,ny,nx,nt);
            if (mode == 3 || mode == 2)
               % w(1:nx*ny*(nt-1),1) = (1 + (scale_const)^2/delta).*w(1:nx*ny*(nt-1),1);
               % w((nt-1)*nx*ny+1:end,1) = ((1/2) + (scale_const)^2/delta) .* w((nt-1)*nx*ny+1:end,1);
                w(1:nx*ny*(nt-1),1) = (beta*(1 + 1/delta)).*w(1:nx*ny*(nt-1),1); % CHANGE HAPPENED HERE.
                w((nt-1)*nx*ny+1:end,1) = (beta*((1/2) + 1/delta)) .* w((nt-1)*nx*ny+1:end,1); % CHANGE HAPPENED HERE.
                x = w + x;
            else
                w((nt-1)*nx*ny+1:end,1) = (1/2) .* w((nt-1)*nx*ny+1:end,1);
                x = beta.*w + x;
            end
        end 

        function x = Apply_DB_Preconditioner(w)
            x = Lvl3_Circulant_Operator(P,w,ny,nx,nt,"inv");
        end
    end

    % ==================================================================================================================== %
    % header + output printing functions: 
    % pl = 1: primal-dual infeasibility and mu is printed at each iteration k
    % pl = 2: primal-dual infeasibility, mu, sigma, and step-lengths are printed at each iteration k
    % -------------------------------------------------------------------------------------------------------------------- %
    function header()
        fprintf(' ');
        fprintf('%4s    ', 'iter');
        fprintf('%8s  ', 'pr feas: Equalities');
        fprintf('%8s  ', 'pr feas: Inequalities');
        fprintf('\n ====  ======================  ======================');
        fprintf('\n');
    end


    function output(it,res_p_eq,res_p_ineq)
        fprintf(' ');
        fprintf('%4d    ', it);
        fprintf('%16.2e  ', res_p_eq);
        fprintf('%16.2e  ', res_p_ineq);
        fprintf('\n');
    end
% ==================================================================================================================== %
% ******************************************************************************************************************** %
% END OF FILE
% ******************************************************************************************************************** %