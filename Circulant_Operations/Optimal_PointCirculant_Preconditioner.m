function [circulant_values] = Optimal_PointCirculant_Preconditioner(A,r)
% ==================================================================================================================== %
% [circulant_values] = Optimal_PointCirculant_Preconditioner(A,r)
% -------------------------------------------------------------------------------------------------------------------- %
% This function takes a dense Toeplitz matrix T as an argument and computes its optimal circulant preconditioner.
% 
% Finally, the function returns a vector that holds all the distinct circulant values
% of the optimal circulant preconditioner.
%
% Input: r = 0 (default), general Toeplitz matrix -> matrix represented as an (2*n-1)-dimensional vector. 
%                                                    This vector should be given in the following form:
%                                                            t = [t_1,...,t_n,t_{-1},...,t_{-n+1}]^T
%        r = 1, symmetric Toeplitz matrix -> matrix represented as an n-dimensional vector.
%        r = 2, general dense matrix -> matrix represented as an n by n matrix.
% Author: Spyridon Pougkakiotis.
% ==================================================================================================================== %

    % ================================================================================================================ %
    % Input check for potential errors
    % ---------------------------------------------------------------------------------------------------------------- %
    if (nargin < 2 || isempty(r))
        r = 0;
    end
    if (size(A,1) ~= size(A,2)) % Either a dimension error, or a Toeplitz vector.
        if (size(A,1) > 1 && size(A,2) > 1)
            error("Matrix is not square.");
        elseif (size(A,2) > 1) 
            A = A'; 
        end
        if (r == 0)
            n = (size(A,1)+1)/2; % The vector is assumed to be of size 2*n-1.
        else
            n = size(A,1);
        end
    elseif (size(A,1) == size(A,2))
         n = size(A,1);
    else
        error("This function accepts only dense Toeplitz matrices, represented as a vector.")
    end
    % ================================================================================================================ %

    % ================================================================================================================ %
    % Building the optimal circulant preconditioner, i.e. the one that solves the following problem:
    %                          C = argmin_{C in C_n} ||C - A||_F,
    % where C_n is the set of Circulant  matrices, of order n.
    % ---------------------------------------------------------------------------------------------------------------- %
    
    % ================================================================================================================ %
    % Dense Toeplitz vector as an input. Forming the optimal circulant preconditioner.
    % ---------------------------------------------------------------------------------------------------------------- %
    circulant_values = zeros(n,1);
    if (r == 0)
        for i = 1:n
            if (i == 1)
                circulant_values(i) = A(i);
            else
                circulant_values(i) = (1/n)*((n-i+1)*A(i) + (i-1)*A(2*n-(i-1)));
            end
        end
    elseif (r == 1)
        for i = 1:n
            if (i == 1)
                circulant_values(i) = A(i);
            else
                circulant_values(i) = (1/n)*((n-i+1)*A(i) + (i-1)*A(n-i+2));
            end
        end
    elseif (r == 2)
        % ======================================================================================================== %
        % General dense matrix as an input. Forming the optimal preconditioner.
        % -------------------------------------------------------------------------------------------------------- %
        circulant_values = zeros(n,1);
        for i = 1:n
            for j = 1:n
                circulant_values(mod(i-j,n)+1) = circulant_values(mod(i-j,n)+1) + A(i,j);
            end
        end
        circulant_values = (1/n) .* circulant_values;
        % ======================================================================================================== %
    else
        error("Wrong structure parameter.\n")
    end
    % ================================================================================================================ %
   
    % ================================================================================================================ %

end
% ******************************************************************************************************************** %
% END OF FILE
% ******************************************************************************************************************** %

