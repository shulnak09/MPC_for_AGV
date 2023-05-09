function [H, f, Aeq, beq, lb, ub] = Prepare_QP_solver(A_d, B_d, d,  x_lb,x_ub, xf_lb,xf_ub,u_lb,u_ub,P, Q, R , N, x0, x_desired)

[n,~] = size(A_d);
[~,m] = size(B_d);

% Define the H matrix:
% Define Q_bar:
Q_bar = zeros(n*N,n*N);
for k = 1:N-1
    Q_bar((k-1)*n+1:k*n,(k-1)*n+1:k*n) = 2*Q;
end
Q_bar(n*(N-1)+1:n*N,n*(N-1)+1:n*N) = 2*P;

% Define R_bar:
R_bar = zeros(m*N, m*N);
for k = 1:N
    R_bar((k-1)*m+1:k*m,(k-1)*m+1:k*m) = 2*R;
end

H = [Q_bar zeros(n*N,m*N)
    zeros(m*N,n*N) R_bar];

% Compute f:
f = zeros(N*(m+n),1);
% f = [-2*Q_bar*reshape(x_desired,[],1); zeros(m*N,1)];
for k = 1:N-1
    f((k-1)*n+1:k*n,1) = -2*Q*x_desired;
end
f(n*(N-1)+1:n*N,1) = -2*P*x_desired;


% Define the equality constraints:
Aeq                  = zeros(N*n,N*(n+m));
Aeq(1:n,1:n)         = eye(n);
Aeq(1:n,n*N+1:n*N+m) = -B_d;
for k=2:N
   Aeq((k-1)*n+1:k*n,(k-1)*n+1:k*n)         = eye(n);
   Aeq((k-1)*n+1:k*n,(k-2)*n+1:(k-1)*n)     = -A_d;
   Aeq((k-1)*n+1:k*n,n*N+(k-1)*m+1:n*N+k*m) = -B_d;
end


beq         = zeros(N*n,1);
beq(1:n, 1) = A_d*x0 + d;
for k =2:N
    beq((k-1)*n+1:k*n,1) = d;
end



% Compute the lower and upper bounds:
x_lb_vec = zeros(n*N,1);
x_ub_vec = zeros(n*N,1);

for k = 1:N-1
    x_lb_vec((k-1)*n+1:n*k,1) = x_lb;
    x_ub_vec((k-1)*n+1:n*k,1) = x_ub;
end

x_lb_vec((N-1)*n+1:n*N,1) = xf_lb;
x_ub_vec((N-1)*n+1:n*N,1) = xf_ub;

u_lb_vec = zeros(m*N,1);
u_ub_vec = zeros(m*N,1);

for k = 1:N
    u_lb_vec((k-1)*m+1:m*k,1) = u_lb;
    u_ub_vec((k-1)*m+1:m*k,1) = u_ub;
end

lb = [x_lb_vec; u_lb_vec];
ub = [x_ub_vec; u_ub_vec];


end

