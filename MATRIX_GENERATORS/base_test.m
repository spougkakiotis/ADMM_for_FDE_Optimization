nx = 8; nt = 8^2; gamma = 1e0;

[M1,M2,M3,A,f1,f2,f3] = FDE_matrices(0,1,0,1,0,1,nx,nx,nt,1.5,1.5,0.5);

Matrix = [M1 sparse(nx^2*nt,nx^2*nt) A';
    sparse(nx^2*nt,nx^2*nt) gamma*M2 M3;
    A M3 sparse(nx^2*nt,nx^2*nt)];
rhs = [f1; f2; zeros(nx^2*nt,1)];

Matrix2 = A'*(M2*A)+10^-4*(M1);

Matrix3  = gamma*A'*(M2*A)+(M1);

sol = Matrix\rhs;
% Plot state solution at last-but-one time-step
figure (1), surf(reshape(sol((nt-2)*nx^2+1:(nt-1)*nx^2),[nx nx]))
% Plot control solution at last-but-one time-step
figure (2), surf(reshape(sol((2*nt-2)*nx^2+1:(2*nt-1)*nx^2),[nx nx]))
% Plot desired state at last-but-one time-step
figure (3), surf(reshape(f1((nt-2)*nx^2+1:(nt-1)*nx^2),[nx nx]))


