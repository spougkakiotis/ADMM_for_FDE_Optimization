function P = NE_preconditioner(A, A_t, circ_const, gamma, n_1, n_2, n_3,mode,M1,delta,scale_const,beta)
% ============================================================================================================================================================ %
% This function takes as an input the data of the FDE optimization problem, 
% and returns the Level-3 circulant preconditioner of the mode equations' matrix, i.e. M + gamma A^T M A. It
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
    
    CA_T_eigs = Lvl3_Circulant_Preconditioner(A_t,n_1,n_2,n_3);
    if (mode == 0 || mode == 1)
        CA_eigs = CA_eigs.*(gamma*circ_const);
        P = Lvl3_Circulant_Mat_Mul(CA_eigs,CA_T_eigs);
        if (mode == 0)
            P = P + (scale_const^2)*circ_const;
        else
            M1_eigs = Lvl3_Circulant_Preconditioner(M1,n_1,n_2,n_3);
            M1_eigs = ((scale_const^2)*(circ_const)).*(M1_eigs);
            P = P + M1_eigs;
        end
    elseif (mode == 2 || mode == 3 || mode == 4)
        if (mode == 2) %inequality only on y
            CA_eigs = CA_eigs.*(1/((scale_const)^2/(beta*gamma*circ_const)+delta));
        elseif (mode == 3)
            CA_eigs = CA_eigs.*(1/((scale_const^2)/(beta*(gamma*circ_const + (1/delta)))+delta/beta));
        else
            CA_eigs = CA_eigs.*(1/((scale_const^2)/(beta*(gamma*circ_const + (1/delta)))+delta/beta)); % Change happened here.
        end
        P = Lvl3_Circulant_Mat_Mul(CA_eigs,CA_T_eigs);
        if (mode ~= 4)                                       % Change happened here.
           % P = P + circ_const + ((scale_const)^2/delta);
            P = P + beta*(circ_const + (1/delta));
        else
            P = P + beta*circ_const;
        end
    end
end

