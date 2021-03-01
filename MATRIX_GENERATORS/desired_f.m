function v = desired_f(t,x,y)

%   v = f(t,x,y);
%   input
%          t            time variable
%          x            x spatial variable
%          y            y spatial variable
%   output
%          v            value of the desired state of the FDE at the point
%                       (x,y) at time t
%
      %  v = sin(pi*x)*sin(pi*y);
        v = 10 * cos(10 * x) * sin(x * y)*(1-exp(-5*t));
      % v = sin(pi*x)*sin(pi*y) + 5*sin(2*pi*x)*sin(3*pi*y) + 3*cos((pi/2)*x)*sin(pi*x);
end