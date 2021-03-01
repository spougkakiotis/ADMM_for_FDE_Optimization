function [product] = PointCirculantOperator(c,x)
% ==================================================================================================================== %
% [product] = PointCirculantOperator(c,x): 
% -------------------------------------------------------------------------------------------------------------------- %
% This function, takes as an input the first column of a circulant matrix, say c, as well as a vector, x.
% Then, using fft and ifft, it efficiently computes the matrix vector product: C*x, where C is the 
% circulant matrix produced by vector c. Then, the result is returned. Overall complexity: O(nlogn)
% More specifically, a circulant matrix C (of dimension n by n), can be written as:
%                         C = F_n^(-1) * diag(F_n*c) * F_n
% Hence, we can compute a matrix vector product C*x as:
%                       product =  F_n^(-1)(F_n(c) .* F_n(x))
% in O(nlogn) operations.
% ==================================================================================================================== %

% ==================================================================================================================== %
% COMMENTS:
% -------------------------------------------------------------------------------------------------------------------- %
% One potential improvement for this code, woudl be to use the function fftw(), which optimizes fft for the specific
% size of the input data, and their structure. However, in order for this to be efficient, fft() must be called several
% times with the same or similar input data. See "help fftw".
% ==================================================================================================================== %


    % ================================================================================================================ %
    % Input check for potential errors
    % ---------------------------------------------------------------------------------------------------------------- %
    if (nargin < 1 || isempty(c))
        error("Not enough input arguments");
    elseif (nargin < 2 || isempty(x))
        error("Not enough input arguments");
    end
    if (size(c,2) > 1) 
        c = c';
    end
    if (size(x,2) > 1)
        x = x'; 
    end
    if (size(c,2) > 1)
        error("Error in data dimension. 1st arg. must be a vector, containing the first column of a circulant matrix.");
    end
    if (size(x,2) > 1)
        error("Error in data dimension. 2nd arg. must be a vector.");
    end
    if (size(x,1) ~= size(c,1))
        error("Incompatible dimensions.");
    end
    % ================================================================================================================ %

    % ================================================================================================================ %
    % Matrix-vector multiplication using fft, ifft
    % ---------------------------------------------------------------------------------------------------------------- %
    lambda = fft(c);       % Compute the eigenvalues of C. 
    temp = fft(x);         % Compute the Fourier transform of x.
    temp = lambda .* temp; 
    product = ifft(temp);  % Finaly, produce the result.
    % ================================================================================================================ %
end
% ******************************************************************************************************************** %
% END OF FILE
% ******************************************************************************************************************** %