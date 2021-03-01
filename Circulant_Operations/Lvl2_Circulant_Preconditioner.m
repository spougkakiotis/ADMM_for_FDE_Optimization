function [mat_eigenvals] = Lvl2_Circulant_Preconditioner(A,n_1,n_2,prec_choice)
% ============================================================================================================================================================ %
% Building the level-2 circulant preconditioner
% This function assumes that A is Doubly-symmetric Doubly-Block Toeplitz (i.e. size(A,1) = n_1, size(A,2) = n_2).
% It returns a matrix containing the eigenvalues of the level-2 Circulant matrix (since the decomposition is 
% constant and hence needs not be computed for different matrices), i.e. 
% it returns a matrix D such that:
%            C_l2 = tilde_F* * D * tilde_F, 
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
    % ================================================================================================================ %

   
    % ======================================================================================================== %
    % Step 1: Apply the Optimal_PointCircunalnt_Preconditioner to each of the 2*m-1 blocks.
    % -------------------------------------------------------------------------------------------------------- %
    for i = 1:n_2
        if (prec_choice == "optimal")
            A(:,i) = Optimal_PointCirculant_Preconditioner(A(:,i),1);
        elseif (prec_choice == "strang")
            A(:,i) = Strang_PointCirculant_Preconditioner(A(:,i),1);
        end
    end
    A = fft(A,[],1);
    % ======================================================================================================== %
        
    % ======================================================================================================== %
    % Step 2: Compute the matrices D_k such that: (D_k)_{ij} = (lambda_{ij})_{kk}, where lambda are the
    %         eigenvalues of each of the m circulant blocks. (k here is k = 1,...,n).
    % -------------------------------------------------------------------------------------------------------- %
    A = A';
    % ======================================================================================================== %
        
    % ======================================================================================================== %
    % Step 3: Compute the circulant approximation of each of the dense block matrices D_k and its eigenvalues.
    % -------------------------------------------------------------------------------------------------------- %
    mat_eigenvals = zeros(n_2,n_1);
    for i = 1:n_1
        if (prec_choice == "optimal")
            mat_eigenvals(:,i) = Optimal_PointCirculant_Preconditioner(A(:,i),1);
        elseif (prec_choice == "strang")
            mat_eigenvals(:,i) = Strang_PointCirculant_Preconditioner(A(:,i),1);
        end
    end
    mat_eigenvals = fft(mat_eigenvals,[],1);
    mat_eigenvals = reshape(mat_eigenvals,n_1*n_2,1); 
    % ======================================================================================================== %
end
% ******************************************************************************************************************** %
% END OF FILE
% ******************************************************************************************************************** %
