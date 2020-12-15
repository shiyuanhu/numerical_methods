function [p_fine, p_nodes] = cubic_spline(f,t_nodes,t_fine)
% compute the natural cubic spline interpolation
% inputs:
%       f: anonymous function to be interpolated
%       t_nodes: node points where the value and the derivative of the
%                 function is known
%       t_fine: a finer grid evaluate the interpolated polynomial.
% 
% Inside iterval i, the cubic polynomial is given by
% S_i(x) = p0 + p1*(x-x_i) + p2*(x-x_i)^2 + p3*(x-x_i)^3,
% where x_i < x < x_{i+1}
%
n = length(t_nodes); % number of node points
h = t_nodes(2:end)-t_nodes(1:end-1);
y = f(t_nodes);
%
% form the sparse matrix
S = matrix_sigma(h,n);
%
% form the right hand side given by Eq. (11.7) in S&M
df = (y(2:n) - y(1:n-1))./h;
rhs = 6*(df(2:end)-df(1:end-1));
%
% solve sigma
sigma = (S\rhs); sigma = [0; sigma; 0];
%
% form coefficients of the polynomials
p0 = y;
p1 = df - h.*(2*sigma(1:end-1) + sigma(2:end))/6;
p2 = sigma./2;
p3 = (sigma(2:end)-sigma(1:end-1))./(6*h);

% form piecewise polynomials
p_fine = interp(t_fine,p0,p1,p2,p3);
p_nodes = interp(t_nodes, p0,p1,p2,p3);

    function results = interp(t,p0,p1,p2,p3)
        % form the cubic piecewise polynomial using Eq. (11.5) in S&M.
        results = zeros(length(t),1); kk = 1;
        for i = 1:n-1
            % find the fine grid for the current interval
            if i ~= n-1
                t_inter = t(t>=t_nodes(i) & t<t_nodes(i+1));
            else
                t_inter = t(t>=t_nodes(i) & t<=t_nodes(i+1));
            end
            nt = length(t_inter);
            % form the cubic polynomial for the current interval
            results(kk:kk+nt-1) = p0(i)+p1(i)*(t_inter-t_nodes(i))+p2(i)*(t_inter-t_nodes(i)).^2+...
                       p3(i)*(t_inter-t_nodes(i)).^3;
            kk = kk+nt;
        end
    end
    
    function S = matrix_sigma(h,n)
        % form the sparse matrix given by Eq. (11.7) in S&M
        onDiags = 2*(h(2:end)+h(1:end-1));
        lower = h(1:end-1); upper = h(2:end);
        S = spdiags([lower onDiags upper], -1:1, n-2, n-2);
    end
end