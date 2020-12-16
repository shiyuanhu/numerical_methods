function x = conju_grad(A,b,max_iter,error_thres)
% conjugate gradient iteration solving the linear system Ax=b. 
% max_iter: maximum interations; error_thres: threshold residual
% rp: r_{n-1}; r: r_{n}
x = 0; rp = b; p = rp; 
n = 0;
residue = inf;
while residue>error_thres
    if n>max_iter % terminate if maximum iterations is reached
        disp('Maximum number of iteration is reached!')
        break
    end
    Ap = A*p;
    alpha = rp'*rp/(p'*Ap);
    if alpha<0 % A much be symmetric positive definite
        error('A must be positive definite!');
    end
    x = x+alpha*p; % update x_n
    r = rp-alpha*Ap; % update residual vector
    residue = norm(r); % residual is the norm
    beta = r'*r/(rp'*rp);
    rp = r; % update r_{n-1}
    p = r+beta*p;
    n = n+1;
    disp(n)
end