function [p_fine, p_nodes] = hermite_inter(f,fderiv,t_nodes,t_fine)
% Hemite interpolation of a function at t_nodes. 
% inputs:
%        f: anonymous function to be interpolated
%        fderiv: derivative of the function f
%        t_nodes: node points where the value and the derivative of the
%                 function is known
%        t_fine: a finer grid evaluate the interpolated polynomial.
%
p_fine = 0; p_nodes = 0;
for j = 1:length(t_nodes)
    p_fine = p_fine+Hk(j,t_fine)*f(t_nodes(j))+Kk(j,t_fine)*fderiv(t_nodes(j));
    p_nodes = p_nodes+Hk(j,t_nodes)*f(t_nodes(j))+Kk(j,t_nodes)*fderiv(t_nodes(j));
end
    
    function coeff = Lk(k,t)
        % evaluate Lk given by Eq. (6.4) in S&M
        % input:
        %       k: the node number
        %       t: evalute Lk at t
        coeff = 1;
        for i=1:length(t_nodes)
            if i ~= k
                coeff = coeff.*((t-t_nodes(i))./(t_nodes(k)-t_nodes(i)));
            end
        end
    end
    
    function coeff = Lk_deriv(k,t_val)
        % evaluate Lk' using symbolic computation
        % or one can compute Lk' analytically, which should be much faster.
        % input: 
        %       k: the node number
        %       t_val: evalute Lk' at t_val
        t = [];
        syms Lt(t)
        Lt(t) = @(t)Lk(k,t);
        df = diff(Lt);
        coeff = double(df(t_val));
    end

    function coeff = Hk(k,t)
        % evaluate Hk
        coeff = Lk(k,t).^2.0.*(1-2*Lk_deriv(k,t_nodes(k))*(t-t_nodes(k)));
    end

    function coeff = Kk(k,t)
        % evaluate Kk
        coeff = Lk(k,t).^2.0.*(t-t_nodes(k));
    end
    
end