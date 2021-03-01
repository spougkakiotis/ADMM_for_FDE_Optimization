function g_alpha_vector = g_alpha_eval(alpha,n)

%   g_alpha_vector = g_alpha_vector(alpha, n);
%   input
%          alpha              exponent of the fractional derivative
%          n                  number of points on the discretized interval
%   output
%          g_alpha_vector     vector containing the coefficient of the
%                             discretized fractional derivative
%


    g_alpha_vector=zeros(n,1);
    
    g_alpha_vector(1,1)=1;
    
    for i = 2:n
        g_alpha_vector(i,1)=(1-(alpha+1)/(i-1))*g_alpha_vector(i-1,1);
    end
    
end