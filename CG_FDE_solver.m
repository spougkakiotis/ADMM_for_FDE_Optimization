function [y,u,p,flag,res,iter] = CG_FDE_solver(A,A_t,g,y_hat,gamma,ny,nx,nt,P,obj_norm,M1,scale_const)
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
    maxit = 350;     % maximum number of CG iterations.
    tol = 10^(-8);  % error tolerance for CG termination.
    lvls = [ny;nx;nt];
    A = reshape(A,ny*nx*(2*nt-1),1);
    A_t = reshape(A_t,ny*nx*(2*nt-1),1);
    Acmat = Multilevel_Circulant_Extrapolation(A, lvls, 3);
    A_tcmat = Multilevel_Circulant_Extrapolation(A_t, lvls, 3);
    if (obj_norm == 1)
        M1 = reshape(M1,ny*nx*(2*nt-1),1);
        M1_cmat = Multilevel_Circulant_Extrapolation(M1, lvls, 3);
    end
    function x = NE_multiplier(w)
        x = Lvl3_Toeplitz_Operator(Acmat,w,ny,nx,nt);
        x((nt-1)*nx*ny+1:end,1) = (1/2) .* x((nt-1)*nx*ny+1:end,1); % Apply operator M2.
        x = Lvl3_Toeplitz_Operator(A_tcmat,x,ny,nx,nt);
        if (obj_norm == 0)                                          % Apply operator M1.
            w((nt-1)*nx*ny+1:end,1) = (1/2) .* w((nt-1)*nx*ny+1:end,1);
        elseif (obj_norm == 1)
            w = Lvl3_Toeplitz_Operator(M1_cmat,w,ny,nx,nt);
            w((nt-1)*nx*ny+1:end,1) = (1/2) .* w((nt-1)*nx*ny+1:end,1);
        end
        x = (scale_const^2).*w + gamma.*x;
    end 

    function x = Apply_DB_Preconditioner(w)
        x = Lvl3_Circulant_Operator(P,w,ny,nx,nt,"inv");
    end
    tmp = (scale_const*gamma).*g;
    tmp((nt-1)*nx*ny+1:end,1) = (1/2) .* tmp((nt-1)*nx*ny+1:end,1);
    if (obj_norm == 1) % Then compute M_1*x only in the H-1 obj_norm case.
        y_hat = Lvl3_Toeplitz_Operator(M1_cmat,y_hat,ny,nx,nt);
        y_hat((nt-1)*nx*ny+1:end,1) = (1/2).* y_hat((nt-1)*nx*ny+1:end,1);
    end
    rhs = Lvl3_Toeplitz_Operator(A_tcmat,tmp,ny,nx,nt) + (scale_const^2).*y_hat;
    
    [y,flag,res,iter] = pcg(@(w) NE_multiplier(w), rhs, tol, maxit,@(w) Apply_DB_Preconditioner(w));
    u = g - (1/scale_const).*Lvl3_Toeplitz_Operator(Acmat,y,ny,nx,nt);
    tmp2 = u; 
    tmp2((nt-1)*nx*ny+1:end,1) = (1/2) .* tmp2((nt-1)*nx*ny+1:end,1);
    p = -(gamma/scale_const).*tmp2;
    
   

end

