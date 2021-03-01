function [M1, M2, M3, A, f1, f2, f3,C_alpha,L_beta1] = FDE_matrices(a,b,c,d,t0,T,nx,ny,nt,betax,betay,alpha)

%   [M1, M2, M3, A, f1, f2, f3] = FDE_matrices(a, b, c, d, t0, T, nx, ny, nt, betax, betay, alpha);
%   input
%          a, b, c, d            [a,b]x[c,d] domain on which is defined the FDE
%          t0,T                  [t0,T] time inteval on which is defined the FDE
%          nx                    number of points on the (interior) grid of the x-axis
%          ny                    number of points on the (interior) grid of the y-axis
%          nt                    number of sub-intervals on the time interval
%          betax                 exponent of the fractional derivative with
%                                respect to the x variable
%          betay                 exponent of the fractional derivative with
%                                respect to the y variable
%          alpha                 exponent of the fractional derivative with
%                                respect to the time variable

%   output
%          M1, M2, M3, A         matrices of the KKT system
%          f1, f2, f3            right hand sides of the KKT system
%


    [Caputo_der_vector, RL_der_vector_x, RL_der_vector_y, f1, f2, f3] = FDE_discretize(a,b,c,d,t0,T,nx,ny,nt,betax,betay,alpha,0);
    
    
    
    C_alpha = zeros(nt, nt);
    
    L_beta1 = zeros(nx, nx);
    
    L_beta2 = zeros(ny, ny);
    
    
    
    Caputo_der_vector_help=[Caputo_der_vector(1), zeros(1, nt-1)];
    
    C_alpha = toeplitz(Caputo_der_vector, Caputo_der_vector_help)
    
    
    
    L_beta1 = toeplitz(RL_der_vector_x)
    
    
    
    L_beta2 = toeplitz(RL_der_vector_y);
    
    cond(L_beta1)
    
    L = kron(L_beta1, eye(ny)) + kron(eye(nx), L_beta2);
    
    
    
    A = kron(C_alpha, eye(nx*ny)) - kron(eye(nt), L);
    
    
    transf = eye(nt);
    
    transf(nt,nt) = 0.5;
    
    M1 = kron(transf, eye(nx*ny));
    
    M2 = M1;
    
    
    M3 = eye(nx*ny*nt);

end