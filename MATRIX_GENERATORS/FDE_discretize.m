function [Caputo_der_vector, RL_der_vector_x, RL_der_vector_y, f1, f3, M_1] = FDE_discretize(a,b,c,d,t0,T,nx,ny,nt,betax,betay,alpha,obj_norm)

%   [Caputo_der_vector, RL_der_vector_x, RL_der_vector_y, f1, f2, f3] = FDE_discretize(a, b, c, d, t0, T, nx, ny, nt, betax, betay, alpha);
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
%          obj_norm                  indicates which obj_norm is used. In particular, obj_norm = 0 assumes L-2 obj_norm, obj_norm = 1
%                                assumes H-1 obj_norm (only on the state variables).
%   output
%          Caputo_der_vector     vector containing the coefficient of the
%                                discretized fractional (Caputo) derivative
%          RL_der_vector_x       vector containing the coefficient of the
%                                discretized fractional (Riemann-Liouville)
%                                derivative with respect to the x variable
%          RL_der_vector_y       vector containing the coefficient of the
%                                discretized fractional (Riemann-Liouville)
%                                derivative with respect to the y variable
%          f1                    vector containing the desired state of the
%                                FDE
%          f2                    vector of all zero, required for
%                                constructing the KKT system
%          f3                    vector containing the initial condition
%                                and its effects for time t>t0
%
%          

    % time step lenght
    deltat=(T-t0)/nt;
    
    % lenght of the intervals on the x and y axis
    hx=(b-a)/(nx+1);
    
    hy=(d-c)/(ny+1);
    
    
    % dimension of the spatial grid
    dim = nx*ny;
    

    % evaluation of the coefficients for the discretized Caputo fractional
    % derivative
    Caputo_der_vector_support = Caputo_derivative(alpha,nt+1);
    
    Caputo_der_vector_support = deltat^(-alpha) .* Caputo_der_vector_support;
    
    Caputo_der_vector = Caputo_der_vector_support(1:nt);
    
    
    % evaluation of the coefficients for the discretized Riemann-Liouville
    % fractional derivative with respect to the x variable
    RL_der_vector_x = RL_derivative_vector(betax,nx);
    
    RL_der_vector_x = hx^(-betax) .* RL_der_vector_x;
    
    %RL_der_vector_x = RL_der_vector_x_support(1:nx);
    
    
    % evaluation of the coefficients for the discretized Riemann-Liouville
    % fractional derivative with respect to the y variable
    RL_der_vector_y = RL_derivative_vector(betay,ny);
    
    RL_der_vector_y = hy^(-betay) .* RL_der_vector_y;
    
   % RL_der_vector_y = RL_der_vector_y_support(1:ny);
    
    
    % evaluation of the desired state
    if (obj_norm == 0)
        y_hat = desired_state_case_1(a, b, c, d, t0, T, nx, ny, nt, hx, hy, deltat);
        M_1 = [];
    elseif (obj_norm == 1)
        [y_hat, M_1] = desired_state_case_2(a, b, c, d, t0, T, nx, ny, nt, hx, hy, deltat);
    end
    
    %evaluation of f1
    f1 = y_hat; % Only because it is zero. Otherwise, for the optimality conditions we must multiply by M_1.
    
    % evaluation of f2
    %f2 = zeros(nx*ny*nt,1);
    
    
    % evaluation of the initial condition y0
    y0 = initial_condition(a, b, c, d, t0, nx, ny);
    
    
    % evaluation of f3
    f3 = zeros(nx*ny*nt,1);
    
    for t=1:nt
        f3((t-1)*dim+1:t*dim) = f3((t-1)*dim+1:t*dim) - Caputo_der_vector_support(t+1).*y0;
    end
    
    
    
end