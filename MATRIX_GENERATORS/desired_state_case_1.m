function y_hat = desired_state_case_1(a, b, c, d, t0, T, nx, ny, nt, hx, hy, deltat)

%   y_hat = desired_state(a, b, c, d, t0, T, nx, ny, nt, hx, hy, deltat);
%   input
%          a, b, c, d       [a,b]x[c,d] domain on which is defined the FDE
%          t0, T            [t0,T] time interval on which is defined the FDE
%          nx               number of points on the (interior) grid of the x-axis
%          ny               number of points on the (interior) grid of the y-axis
%          nt               nt+1 is the number of points on the discretized
%                           time interval
%          hx               lenght of each interval on the x-axis
%          hy               lenght of each interval on the y-axis
%          deltat           time step lenght
%   output
%          y_hat            vector containing the desired state of the FDE
%


    y_hat = zeros(nx*ny*nt,1);
    
    
    x = linspace(a,b,nx+2);
    
    y = linspace(c,d,ny+2);
    
    t = linspace(t0,T,nt+1);

    
    x1 = x(2:nx+1);
    
    y1 = y(2:ny+1);
    
    t1 = t(2:nt+1);
    
    
    for k= 1 : nt-1
        for i= 1 : nx
            for j= 1 : ny
                y_hat((k-1)*(nx*ny)+(i-1)*ny+j) = desired_f(t1(k),x1(i),y1(j));
            end
        end 
    end
    
    
	for i= 1 : nx
        for j= 1 : ny
            y_hat((nt-1)*(nx*ny)+(i-1)*ny+j) = 0.5 * desired_f(t1(nt),x1(i),y1(j));
        end
	end 

    
end