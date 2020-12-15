function R = romberg_int(f,a,b,m,nmax,eps)
% Form Romberg integration
% inputs:
%        f: the function handle to be integrated
%        a, b: lower and upper integration limit
%        m: the number of intervals
%        nmax: maximum number of rows in the Romberg table
%        eps: desired accuracy
% returns:
%        R: Romberg table
%           the most accurate result R_{n,n}
%
% Written by Shiyuan Hu <shiyuan.hu@nyu.edu>, Dec. 1 2020
%
% set default accuracy
if nargin<6
    eps = 1e-10;
end
% initialize error and Romberg table R
error = inf;
R = zeros(nmax, nmax);
%
% initialize grid points
x = linspace(a,b,m+1); dx = (b-a)/m;
%
% compute R_{1,1}
fx = f(x);
R(1,1) = trapez(fx,dx);

% iterate through successive doubling until maximum number of
% doubling or desired accuracy is reached
for i = 2:nmax
    if error<eps
        break
    end
    [x,x1,dx] = doubling(x,dx); % double grid points
    % reuse R_{i-1,1} to compute R_{i,1}, and only 
    % evaluate f at new grids
    R(i,1) = 1/2*R(i-1,1)+dx*sum(f(x1)); 
    % Richardson extrapolation
    for k = 2:i
        R(i,k) = R(i,k-1)+(R(i,k-1)-R(i-1,k-1))/(4^(k-1)-1);
    end
    % estimate error
    error = abs(R(i,i)-R(i-1,i-1));
    if i == nmax
        disp('Maximum number of doubling is reached!')
    end
end
% remove zero rows and columns
R(all(~R,2),:) = [];
R(:,all(~R,1)) = [];

    function u = trapez(ux,dx)
        % form trapezoidal integration 
        u = dx*(sum(ux)-0.5*ux(1)-0.5*ux(end));
    end

    function [x,x1,dx] = doubling(x,dx)
        % doubling the grid
        % input:
        %       x: grid points before doubling
        %       dx: grid size before doubling
        % return:
        %       x: grid points after doubling
        %       x1: extra grid points
        %       dx: grid size after doubling      
        x1 = (x(1:end-1)+x(2:end))./2;
        x = sort([x x1]); dx = dx/2.;
    end
end