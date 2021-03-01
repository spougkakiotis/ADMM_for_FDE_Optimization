function [P,M1_eigs] = AS_preconditioner(A, A_t, circ_const, gamma, n_1, n_2, n_3,obj_norm,M1,scale_const)
% ============================================================================================================================================================ %
% This function takes as an input the data of the FDE optimization problem, 
% and returns the Level-3 circulant preconditioner of the obj_normal equations' matrix, i.e. M + gamma A^T M A. It
% assumes that matrix A is Level-3 Toeplitz, represented as such (i.e. size(A,1) = n,
% size(A,2) = (2*nt -1)).
% ============================================================================================================================================================ %
    if (nargin < 3 || isempty(circ_const))
        error("Not enough input arguments.")
    end
    if (nargin < 4 || isempty(gamma))
        error("Not enough input arguments.")
    end
    if (nargin < 5 || isempty(n_1))
        error("Not enough input arguments.")
    end
    if (nargin < 6 || isempty(n_2))
        error("Not enough input arguments.")
    end
    if (nargin < 7 || isempty(n_3))
        error("Not enough input arguments.")
    end
    CA_eigs = Lvl3_Circulant_Preconditioner(A,n_1,n_2,n_3);  
    if (obj_norm == 0)
        CA_eigs = CA_eigs.*(1/circ_const);
        M1_eigs = [];
    elseif (obj_norm == 1)
        M1_eigs = Lvl3_Circulant_Preconditioner(M1,n_1,n_2,n_3);
        M1_eigs = circ_const.*M1_eigs;
        CA_eigs = (CA_eigs)./M1_eigs;
    end
    CA_T_eigs = Lvl3_Circulant_Preconditioner(A_t,n_1,n_2,n_3);
    P = Lvl3_Circulant_Mat_Mul(CA_eigs,CA_T_eigs);
    P = P + (scale_const)^2/(gamma*circ_const);
end

