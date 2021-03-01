function [y_hat, M_1] = desired_state_case_2(a, b, c, d, t0, T, nx, ny, nt, hx, hy, deltat)

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
%          y_hat            vector containing the desired state of the FDE,
%                           modified in such a way that the norm in the
%                           functional is in H^1
%          K1               vector containing the coefficient of the
%                           matrix defined in the functional (first blocks)
%          K2               vector containing the coefficient of the
%                           matrix defined in the functional (last block)
%


    y_hat = zeros(nx*ny*nt,1);
    
    
    x = linspace(a,b,nx+2);
    
    y = linspace(c,d,ny+2);
    
    t = linspace(t0,T,nt+1);

    
    x1 = x(2:nx+1);
    
    y1 = y(2:ny+1);
    
    t1 = t(2:nt+1);
    
    
    for k= 1 : nt%-1
        for i= 1 : nx
            for j= 1 : ny
                y_hat((k-1)*(nx*ny)+(i-1)*ny+j) = desired_f(t1(k),x1(i),y1(j));
            end
        end 
    end
    
    K_x = zeros(nx,1);
    K_x(1,1) = 2*(hx^(-2));
    K_x(2,1) = (-1)*(hx^(-2));
    K_y = zeros(ny,1);
    K_y(1,1) = 2*(hy^(-2));
    K_y(2,1) = (-1)*(hy^(-2));

    K = zeros(ny*nx,1); % doubly-symmetric, doubly-block Toeplitz
    K(1:ny,1) = K_y;
    K(1,1) = K(1,1) + K_x(1,1);
    K(ny+1,1) = K_x(2,1);
    % This is M_1, without the (1/2) in the end. This should be incorporated after the matrix-vector product.
    M_1 = zeros(nx*ny,2*nt-1);
    M_1(:,1) = K;
    M_1(1,1) = M_1(1,1) + 1; %treat this as a level-3 Toeplitz, calculate M_1*x and then multiply the last nx*ny elements of the result by (1/2)
end