function w = Lvl3_Circulant_Operator(C,x,n_1,n_2,n_3,operation)
% ============================================================================================================================================================ %
% w = DoublyBlock_Circulant_Operator(C,x,m)
% ------------------------------------------------------------------------------------------------------------------------------------------------------------ %
% This function takes as an input a level-3 circulant matrix C, represented by its eigenvalues, as well as a vector x, 
% and returns the result of the multiplication (inversion) C*x(C^{-1}*x), in O(n_1 n_2 n_3 log (n_1 n_2 n_3)) operations, where, 
% n_1 - the level 1 size and n_2 - the level 2 size, and so on. 
% The user is able to:
%                     operation = "mul" (default), multiply C by the vector,
%                     operation = "inv", multiply C^{-1} by the vector.
% ____________________________________________________________________________________________________________________________________________________________ %
    if (nargin < 3 || isempty(n_1))
        error("Not enough input arguments.")
    end
    if (nargin < 4 || isempty(n_2))
        error("Not enough input arguments.")
    end
    if (nargin < 5 || isempty(n_3))
        error("Not enough input arguments.")
    end
    if (nargin < 6 || isempty(operation))
        operation = "mul";
    end
    if (size(x,1) ~= n_1*n_2*n_3)
        error("Incorrect input dimensions.")
    end
    
    % ============================================================================================ %
    % Apply the level-3 fourier-like operator, to each of the n_3 blocks of x
    % -------------------------------------------------------------------------------------------- %
    w = zeros(n_1*n_2,n_3);
    for i = 1:n_3
        x_tmp = reshape(x((i-1)*(n_1*n_2)+1:i*(n_1*n_2),1),n_1,n_2);
        x_tmp = fft(x_tmp,[],1);
        x_tmp = x_tmp';
        x_tmp = fft(x_tmp,[],1);
        w(:,i) = reshape(x_tmp,n_1*n_2,1);
    end
    w = w';
    w = fft(w,[],1);
    w = reshape(w,n_1*n_2*n_3,1);
    % ___________________________________________________________________________________________ %
    
    
    % =========================================================================================== %
    % Apply the operation to the eigenvalues of the matrix.
    % ------------------------------------------------------------------------------------------- %
    if (operation == "mul")
        w = C.*w;
    elseif (operation == "inv")
        w = w./C;
    end
    % __________________________________________________________________________________________ %
    
    % ========================================================================================== %
    % Apply the inverse level-3 Fourier-like operator to get the final result.
    % ------------------------------------------------------------------------------------------ %
    w = reshape(w,n_3,n_1*n_2);
    w = ifft(w,[],1);
    w = w';
    for i = 1:n_3
        x_tmp = reshape(w(:,i),n_2,n_1);
        x_tmp = ifft(x_tmp,[],1);
        x_tmp = x_tmp';
        x_tmp = ifft(x_tmp,[],1);
        w(:,i) = reshape(x_tmp,n_1*n_2,1);
    end
    w = real(reshape(w,n_1*n_2*n_3,1));
    % __________________________________________________________________________________________ %
end

