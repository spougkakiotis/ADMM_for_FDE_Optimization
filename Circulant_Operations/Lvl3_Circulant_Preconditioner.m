function [mat_eigenvals] = Lvl3_Circulant_Preconditioner(A,n_1,n_2,n_3,prec_choice)
% ============================================================================================================================================================ %
% Building the level-3 circulant preconditioner.
% This function assumes that A is level-3 Toeplitz with doubly-block doubly-symmetric blocks (i.e. size(A,1) = n_1, size(A,2) = n_2, size(A,3) = 2*n-3-1).
% It returns a matrix containing the eigenvalues of the level-3 Circulant matrix (since the decomposition is 
% constant and hence needs not be computed for different matrices), i.e. 
% it returns a matrix D such that:
%            C_l3 = tilde_F* * D * tilde_F, 
% where tilde_F is an appropriate matrix.
% ============================================================================================================================================================ %

    
    % ================================================================================================================ %
    % Input check for potential errors
    % ---------------------------------------------------------------------------------------------------------------- %
    if (nargin < 2 || isempty(n_1)) 
        error("Not enough input arguments.")
    end
    if (nargin < 3 || isempty(n_2))
        error("Not enough input arguments.")
    end
    if (nargin < 4 || isempty(n_3))
        error("Not enough input arguments.")
    end
    if (nargin < 5 || isempty(prec_choice))
        prec_choice = "optimal";
    end
    if (prec_choice ~= "optimal" && prec_choice ~= "strang")
        error("Choice of preconditioner parameter error");
    end
    
    % ================================================================================================================ %
    lvl2_eigen_blocks = zeros(n_2*n_1,2*n_3-1);
    A = reshape(A,n_1,n_2,2*n_3-1);
    for i = 1:(2*n_3-1)
        lvl2_eigen_blocks(:,i) = Lvl2_Circulant_Preconditioner(A(:,:,i),n_1,n_2,prec_choice);
    end
    lvl2_eigen_blocks = lvl2_eigen_blocks';
    mat_eigenvals = zeros(n_3,n_1*n_2);
    for i = 1:(n_1*n_2)
        if (prec_choice == "optimal")
            mat_eigenvals(:,i) = Optimal_PointCirculant_Preconditioner(lvl2_eigen_blocks(:,i),0);
        elseif (prec_choice == "strang")
            mat_eigenvals(:,i) = Strang_PointCirculant_Preconditioner(lvl2_eigen_blocks(:,i),0);
        end
    end
    mat_eigenvals = fft(mat_eigenvals,[],1);
    mat_eigenvals = reshape(mat_eigenvals,n_1*n_2*n_3,1);
end
% ******************************************************************************************************************** %
% END OF FILE
% ******************************************************************************************************************** %
