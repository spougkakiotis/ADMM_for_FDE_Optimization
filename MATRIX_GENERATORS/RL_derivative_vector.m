function RL_der_vector = RL_derivative_vector(beta,n)

%   RL_der_vector = RL_derivative_vector(beta, n);
%   input
%          beta             exponent of the FDE
%          n                number of sub-intervals on the interval of the x-(y-)axis
%   output
%          RL_der_vector    vector containing the coefficients of the
%                           discretized fractional (Riemann-Liouville) derivative
%


    RL_der_vector = zeros(n,1);
    
    RL_der_vector2 = g_alpha_eval(beta,n+1);
    
    RL_der_vector(3:end,1) = RL_der_vector2(4:end,1);
    
    RL_der_vector(2,1) = RL_der_vector2(1,1) + RL_der_vector2(3,1);
    
    RL_der_vector(1,1) = 2 * RL_der_vector2(2,1);
    
    RL_der_vector = -(0.5*(1/cos(beta*pi/2))) .* RL_der_vector; % Here I have included the Riesz potential scaling term.
  %  RL_der_vector = (0.5) .* RL_der_vector; 
end