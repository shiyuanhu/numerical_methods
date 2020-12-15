% Solve a nonlinear ordinary differential equation using Newton's method
% u^{k+1} = u^{k}-J^{-1}*F(u^{k}), where J(u) is the Jacobian matrix and 
% the nonlinear equation is F(u) = 0. In this example, a nonlinear ODE is 
% solved. -u''(t)-e^{u(t)} = 0, for t in [0,1]. The boundary conditions 
% are u(0) = 0, and u(1) = 0. 
%
global npts h
% subdivide the interval
npts = 128; % number of grid points
h = 1./(npts+1); t = 0:h:1;
%
% initialize u
u = 0.5*t(2:end-1).*(t(2:end-1)-1);
%
% set error tolerence and maximum number of iteration
eps = 1e-15; error = inf; 
nmax = 10; count = 0;
%
% Newton's iteration
while error>eps
    if count>nmax
        disp('Maximum iteration reached!');
        break
    end
    J = Jacobian(u); F = Rhs(u);
    delta_u = (J\(-F'))';
    u = u+delta_u;
    error = max(abs(delta_u));
    count = count+1;
end

function J = Jacobian(u)
% form a sparse Jacobian matrix
    global npts h
    ondiag = 2./h^2-exp(u);
    offdiag = -1./h^2;
    D = sparse(1:npts,1:npts,ondiag,npts,npts);
    E = sparse(2:npts,1:npts-1,offdiag*ones(1,npts-1),npts,npts);
    J = D+E+E';
end

function F = Rhs(u)
% form F(u) = 0
    global npts h
    F(1) = -1./h^2*(u(2)-2*u(1))-exp(u(1));
    F(npts) = -1./h^2*(-2*u(npts)+u(npts-1))-exp(u(npts));
    F(2:npts-1) = -1./h^2*(u(3:npts)-2*u(2:npts-1)+u(1:npts-2))-...
                  exp(u(2:npts-1));
end