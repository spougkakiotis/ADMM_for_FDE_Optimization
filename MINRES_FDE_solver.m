function [y,u,p,flag,res,iter] = MINRES_FDE_solver(A,A_t,g,y_hat,gamma,ny,nx,nt,P,circ_const,obj_norm,M1,M1_eigs,scale_const)
% ============================================================================================================================================================ %
% This function takes an FDE-equality constraint optimization problem and
% attempts to solve its optimality conditions directly, using 
% Preconditioned CG. Since such problems possess a symmetric-block doubly-block Toeplitz structure,
% we employ a double-block circulant preconditioner. All the operations substitute 
% matrices by smaller matrices of vectors, due to the structure of the optimization problems in hand. 
% In particular, we exploit diagonal structure, Toeplitz/circulant structure and
% employ special "matrix"-vector multiplication routines.
% ------------------------------------------------------------------------------------------------------------------------------------------------------------ %
% Input: 
%        A, the symmetric-block doubly-block Toeplitz constraint matrix, represented by a matrix of dimensions (nx*ny) by (2*nt-1), containing
%           all nedded information.
%        g, the right hand side of the constraints (nx*ny*nt dimension).
%        y_hat, the desired state.
% Output: 
%        y, u, p, the proposed solution to the problem.
% ____________________________________________________________________________________________________________________________________________________________ %
    maxit = 250;     % maximum number of CG iterations.
    tol = 10^(-6);  % error tolerance for CG termination.
    n = nx*ny*nt;
    lvls = [ny;nx;nt];
    A = reshape(A,ny*nx*(2*nt-1),1);
    A_t = reshape(A_t,ny*nx*(2*nt-1),1);
    Acmat = Multilevel_Circulant_Extrapolation(A, lvls, 3);
    A_tcmat = Multilevel_Circulant_Extrapolation(A_t, lvls, 3);
    
    if (obj_norm == 1)
        M1 = reshape(M1,ny*nx*(2*nt-1),1);
        M1_cmat = Multilevel_Circulant_Extrapolation(M1, lvls, 3);
    end
    
    function x = AS_multiplier(w)
        x = w;
        if (obj_norm == 0)
            x((nt-1)*nx*ny+1:n,1) = (1/2) .* x((nt-1)*nx*ny+1:n,1);                         % Apply operator M1 on x_1.
        else
            x(1:n,1) = Lvl3_Toeplitz_Operator(M1_cmat,w(1:n,1),ny,nx,nt);
            x((nt-1)*nx*ny+1:n,1) = (1/2) .* x((nt-1)*nx*ny+1:n,1);                         % Apply operator M1 on x_1.
        end
        x(1:n,1) = x(1:n,1) + Lvl3_Toeplitz_Operator(A_tcmat,w(2*n+1:3*n,1),ny,nx,nt);   
        x(n+(nt-1)*nx*ny+1:2*n,1) = (1/2) .* x(n+(nt-1)*nx*ny+1:2*n,1);                     % Apply operator M2 on x_2.
        x(n+1:2*n,1) = gamma.*x(n+1:2*n,1) + scale_const.*w(2*n+1:3*n,1);
        x(2*n+1:3*n,1) = Lvl3_Toeplitz_Operator(Acmat,w(1:n,1),ny,nx,nt) + scale_const.*w(n+1:2*n,1);
    end 

    function w = Apply_DB_Preconditioner(w)
        if (obj_norm == 0)
            w(1:n,1) = (1/circ_const).*w(1:n,1);
        else
            w(1:n,1) = Lvl3_Circulant_Operator(M1_eigs,w(1:n,1),ny,nx,nt,"inv");
        end
        w(n+1:2*n,1) = (1/(circ_const*gamma)).*w(n+1:2*n,1);
        w(2*n+1:3*n,1) = Lvl3_Circulant_Operator(P,w(2*n+1:3*n,1),ny,nx,nt,"inv");
    end
    rhs = zeros(3*n,1);
    if (obj_norm == 1) % Then compute M_1*x only in the H-1 obj_norm case.
        y_hat = Lvl3_Toeplitz_Operator(M1_cmat,y_hat,ny,nx,nt);
        y_hat((nt-1)*nx*ny+1:end,1) = (1/2).* y_hat((nt-1)*nx*ny+1:end,1);
    end
    rhs(1:n,1) = y_hat;
    rhs(2*n+1:3*n,1) = scale_const.*g;
    [sol,flag,res,iter] = minres(@(w) AS_multiplier(w), rhs, tol, maxit,@(w) Apply_DB_Preconditioner(w));
    y = sol(1:n,1);
    u = sol(n+1:2*n,1);
    p = sol(2*n+1:3*n,1);
    
end

