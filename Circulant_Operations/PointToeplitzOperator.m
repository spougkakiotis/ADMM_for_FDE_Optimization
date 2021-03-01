function [product] = PointToeplitzOperator(t,x,r)
% ==================================================================================================================== %
% [product] = PointToeplitzOperator(t,x): 
% -------------------------------------------------------------------------------------------------------------------- %
% This function, takes as an input the first column of a Toeplitz matrix, say t, as well as a vector, x.
% Then, using fft, it efficiently computes the matrix vector product: T*x, where T is the 
% Toeplitz matrix produced by vector t. Then, the result is returned. Overall complexity: O(nlogn)
% More specifically, a circulant matrix T (of dimension n by n), can be embedded in a circulant matrix, say C,
% of size 2n by 2n. This means that: Tx = [I 0]*C*[x;0]. The matrix-vector product involving C, can be computed
% in O(2nlog(2n)) operations. Hence, we will call the PointCirculantOperator function 
% (see "help PointCirculantOperator"), to perform the multiplication and then return the result.
% 
% Input: t, a vector representing the Toeplitz matrix,
%        x, a vector with which we want to multiply the matrix,
%        r = 0, (default) assumes that no symmetry is present in the matrix and expects 
%                         the vector t to contain the first column and then the first row of 
%                         the Toeplitz matrix under consideration, that is:
%                                          t = [t_0,t_1,...,t_{n-1},t_{-1},...,t_{-n+1}]^T.
%        r = 1, assumes that the matrix T is symmetric, and expects a vector t, containing only n values.
% ==================================================================================================================== %

    % ================================================================================================================ %
    % Input check for potential errors
    % ---------------------------------------------------------------------------------------------------------------- %
    if (nargin < 1 || isempty(t))
        error("Not enough input arguments");
    elseif (nargin < 2 || isempty(x))
        error("Not enough input arguments");
    end
    if (nargin < 3 || isempty(r))
        r = 0;
    end
    if (size(t,2) > 1) 
        t = t';
    end
    if (size(x,2) > 1)
        x = x'; 
    end
    if (size(t,2) > 1)
        error("Error in data dimension. 1st arg. must be a vector.");
    end
    if (size(x,2) > 1)
        error("Error in data dimension. 2nd arg. must be a vector.");
    end
    n = size(x,1);
    if (size(t,1) ~= 2*n - 1 && r == 0)
        error("Incompatible dimensions.");
    elseif (size(t,1) ~= n && r == 1)
        error("Incompatible dimensions.");
    end
    % ================================================================================================================ %
    if (r == 0)
        % ================================================================================================================ %
        % Embedding the Toeplitz matrix in a circulant one.
        % Create a 2n by 2n circulant matrix, containing all Toeplitz values. Increase the dimensions of vector x also.
        % ---------------------------------------------------------------------------------------------------------------- %
        t_circ = [t(1:n,1); 0; t(2*n-1:-1:n+1,1)]; 
        % ================================================================================================================ %
    elseif (r == 1)
        % ================================================================================================================ %
        % Embedding the Toeplitz matrix in a circulant one.
        % Create a 2n by 2n circulant matrix, containing all Toeplitz values. Increase the dimensions of vector x also.
        % ---------------------------------------------------------------------------------------------------------------- %
        t_circ = [t(1:n,1); 0; t(n:-1:2,1)]; 
        % ================================================================================================================ %
    else
        error("Incorrect structure parameter input.\n")
    end
    temp = [x; zeros(n,1)];
    % Call the method that performs point-circulant matrix-vector products efficiently.
    temp = PointCirculantOperator(t_circ,temp);
    % Keep only the first n entries of the result.
    product = temp(1:n,1);
end
% ******************************************************************************************************************** %
% END OF FILE
% ******************************************************************************************************************** %
