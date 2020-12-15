function I = gauss_quad(f,a,b,n)
% Form gaussian quadrature integration using legnpts from Chebfun
% inputs:
%        f: the function handle to be integrated
%        a, b: lower and upper integration limit
%        n: the number of nodes
% returns:
%        I: the integration value
%
[nodes,weights] = legpts(n,[a,b]);
I = sum(weights.*f(nodes'));
end
