function [y_hat, A, A_t, g, M_1] = Optimization_Problem_Generator(a,b,c,d,t0,T,nx,ny,nt,betax,betay,alpha,obj_norm,scale_const)
% ============================================================================================================================================================ %
% This function takes as an input the specifications of the FDE problem under
% consideration, and generates all the relevant matrices, in the optimization
% format. It takes as an input the following parameters:
% ------------------------------------------------------------------------------------------------------------------------------------------------------------ %
% INPUT:  a,b to determine the interval in dimension 1,
%         c,d to determine the interval in dimension 2,
%         t0, T, the time interval,
%         nx, the number of points in the discretization of the 1st dimension,
%         ny, the number of points in the discretization of the 2nd dimension,
%         nt, the number of points in the discretization of the time dimension,
%         betax, exponent of the fractional derivative with respect to the x variable,
%         betay, exponent of the fractional derivative with respect to the y variable,
%         alpha, exponent of the fractional derivative with respect to the time variable.
% OUTPUT: 
%         M, the matrix representing the Hessian,
%         y_bar, the vector containing the desired state.
%         A, the matrix representing the constraint matrix corresponding to variable y,
%         A_t,
%         g, the right hand side of the constraints,
% ============================================================================================================================================================ %
    [Caputo_der_vector, RL_der_vector_x, RL_der_vector_y, f1, f3, M_1] = FDE_discretize(a,b,c,d,t0,T,nx,ny,nt,betax,betay,alpha,obj_norm);
    L = zeros(nx*ny,1); % Only nx*ny elements needed, since we know it is doubly-symmetric and doubly-Toeplitz.
    L(1:ny,1) = RL_der_vector_y;
    for i = 1:nx
        L((i-1)*ny + 1,1) = L((i-1)*ny + 1,1) + RL_der_vector_x(i,1);
    end
    A = zeros(nx*ny,2*nt-1);
    A_t = zeros(nx*ny,2*nt-1);
    for i = 1:nt
        A(1,i) = Caputo_der_vector(i); 
    end
    A(:,1) = A(:,1) - L;
    A = scale_const .*A;  % This matrix is called B in the paper.
    A_t(:,1) = A(:,1);
    for i = 1:(nt-1)
        A_t(:,i+nt) = A(:,i+1);
    end
    g = f3;
    y_hat = f1;
    
       
end

