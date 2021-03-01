function Caputo_der_vector = Caputo_derivative(alpha,nt)

%   Caputo_der_vector = Caputo_derivative(alpha, nt);
%   input
%          alpha                   exponent of the FDE
%          nt                      number of sub-intervals on the (discretized) time interval
%   output
%          Caputo_der_vector       vector containing the coefficient of the
%                                  discretized fractional (Caputo) derivative
%


    %Caputo_der_vector = zeros(nt,1);
    
    Caputo_der_vector = g_alpha_eval(alpha,nt);

end