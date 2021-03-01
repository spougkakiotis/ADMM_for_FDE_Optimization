% ============================================================================================================================================================ %
% This script contains the set-up parameters of the FDE optimization problem,
% calls functions to construct the discretized optimization problem and then,
% calls a solution method to solve it.
% ------------------------------------------------------------------------------------------------------------------------------------------------------------ %
% The following parameters have to be specified:
% [a,b] -> Interval over which we optimize in dimension x.
% [c,d] -> Interval over which we optimize in dimension y.
% t0    -> Initial time.
% nx    -> Number of points in the discretization along dimension 1.
% ny    -> Number of points in the discretization along dimension 2.
% nt    -> Number of points in the discretization along time dimension.
% betax -> exponent of the fractional derivative with respect to the x variable.
% betay -> exponent of the fractional derivative with respect to the y variable.
% alpha -> exponent of the fractional derivative with respect to the time variable.
% gamma -> regularization parameter.
% ============================================================================================================================================================ %
clc;
clear;
fid = fopen('Experiment_ADMM_norm_results.txt','a+');
% ============================================================================================================================================================ %
% Initializing dimension intervals.
% ------------------------------------------------------------------------------------------------------------------------------------------------------------ %
%Problem # 1 (Easy case)
 %a = -1; b = 1; c = -1; d = 1; t0 = 0; T = 2;   
%Problem # 2 (Hard case)
  a = 0; b = 1; c = 0; d = 1; t0 = 0; T = 1;
% ============================================================================================================================================================ %

% ============================================================================================================================================================ %
% Initializing problem dimensions (number of discretization points).
% ------------------------------------------------------------------------------------------------------------------------------------------------------------ %
l =  7;
nx = 2^l;
ny = 2^l;
nt = ceil(2^(l));
nx = 50;
ny = 50;
nt = ceil(50^(1));
n = nx*ny;
N = nx*ny*nt;
% time step lenght
tau =(T-t0)/nt;   
% lenght of the intervals on the x and y axis
hx=(b-a)/(nx+1);
hy=(d-c)/(ny+1);


% ============================================================================================================================================================ %

% ============================================================================================================================================================ %
% Initializing FDE parameters.
% ------------------------------------------------------------------------------------------------------------------------------------------------------------ %
alpha = 0.7;
betax = 1.3;
betay = 1.3;
gamma = 1e-4;
mode = 3;
beta = 1.618; % Step-length of ADMM. Standard ADMM converges for values beta in (0,1.618)
scale_const =  min(tau^(alpha),hx^(betax)); % Scaling the discretized matrices
if (mode == 2)
    delta = 0.4;
else
    scale_const =  min(tau^(alpha),hx^(betax));
    delta = 10;
end

% ============================================================================================================================================================ %

% ============================================================================================================================================================ %
% Create the optimization problem given the experiment parameters. We get back the following data:
% 1) M = M_1 = M_2, a vector representing the Hessian corresponding to variables y, u (i.e. M in R^(n,1)).
% 2) A, a matrix representing the doubly-block Toeplitz of the constraint matrix (rows = blocks, columns = num of blocks), corresponding to vector y. The 
%       matrix corresponding to vector u is simply an identity matrix.
% 3) y_hat: vector containing the desired state of the FDE.
% 4) g: the right hand side of the constraints.
% ------------------------------------------------------------------------------------------------------------------------------------------------------------ %
[y_hat, A, A_t, g, M1] = Optimization_Problem_Generator(a,b,c,d,t0,T,nx,ny,nt,betax,betay,alpha,0,scale_const);
% ============================================================================================================================================================ %
time = 0; tic;
circ_const = ((n*(nt-1)*1+n*(1/2))/(n*nt)); % Representing matrix M1,M2.
% N = 50^4, opt at: max(y) = 4.83, min(y) = -5.91, max(u) = 491, min(u) = -353
if (mode == 2) % tests: [-6,6], [-5,5], [-4,4], [-2,2] l = 4
    lb = -(1).*(ones(N,1));
    ub = (1).*(ones(N,1));
    u_a = -Inf;
    u_b = Inf;
    y_a = (lb(1));
    y_b = (ub(1));
elseif (mode == 3) % tests: [-150,200] [-4.5,4.5], [-100,150] [-4,4], [-400,400] [-5,5], [-3,3.8] [-100,100] l =4
    lb = -(4).*(ones(2*N,1));
    lb(N+1:end,1) = -350.*ones(N,1);
    ub = (4).*(ones(2*N,1));
    ub(N+1:end,1) = 350.*ones(N,1);
    y_a = (lb(1));
    y_b = (ub(1));
    u_a = (lb(N+1));
    u_b = ub(N+1);
elseif (mode == 4) %[-400,400], [-150,200], [-100,150], [-100,100], [-50,50] l = 4
    lb = -(100).*(ones(N,1));
    ub = (100).*(ones(N,1));
    y_a = -Inf;
    y_b = Inf;
    u_a = (lb(1));
    u_b = (ub(1));
end
P = NE_preconditioner(A, A_t, circ_const, gamma, ny, nx, nt,mode,M1,delta,scale_const,beta);
[y,u,p,w_y,w_u,opt,iter,avg_iter] = Box_Constrained_ADMM(A,A_t,g,y_hat,gamma,ny,nx,nt,P,mode,lb,ub,delta,scale_const,beta);
time = time+toc;
[err,d_inf] = Solution_Statistics(A_t,y,u,p,y_hat,tau,hx,hy,gamma,ny,nx,nt,mode,M1,scale_const,w_y,w_u);
fprintf(fid," Problem Specifics:\n    nx = %d, gamma = %.2e, alpha = %.2f, betax = %.2f.\n Bounds: \n    y_a = %.2f, y_b = %.2f, u_a = %.2f, u_b = %.2f.\n Convergence Statistics:\n    Time = %.2f, ||(y-y_hat)|| = %e, Dual Inf. = %e, CG Iter. = %d, ADMM Iter. = %d, delta = %.2f.\n \n",...
        nx,gamma,alpha,betax,y_a,y_b,u_a,u_b,time,err,d_inf,avg_iter,iter,delta);
fclose(fid);
fprintf(" Problem Specifics:\n     nx = %d, gamma = %.2e, alpha = %.2f, betax = %.2f.\n Bounds: \n    y_a = %.2f, y_b = %.2f, u_a = %.2f, u_b = %.2f.\n Convergence Statistics:\n    Time = %.2f, ||(y-y_hat)|| = %e, Dual Inf. = %e, CG Iter. = %d, ADMM Iter. = %d, delta = %.2f.\n \n",...
        nx,gamma,alpha,betax,y_a,y_b,u_a,u_b,time,err,d_inf,avg_iter,iter,delta);