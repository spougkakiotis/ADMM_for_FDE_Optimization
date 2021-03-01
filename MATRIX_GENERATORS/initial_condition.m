function y0 = initial_condition(a, b, c, d, t0, nx, ny)

%   y0 = initial_condition(a, b, c, d, t0, nx, ny);
%   input
%          a, b, c, d       [a,b]x[c,d] domain on which is defined the FDE
%          t0               time t0 for the FDE
%          nx               number of points on the (interior) grid of the x-axis
%          ny               number of points on the (interior) grid of the y-axis
%   output
%          y0               vector containing the initial condition of the
%                           FDE
%


	y0 = zeros(nx*ny,1);
    
end