function [p_fine,p_nodes] = lagrange_inter(f,t_nodes,t_fine)
% lagrange interpolation
% input: 
%    f: anonymous function to be interpolated: set this to
%               @(x)1./(1+25*x.^2) to see Runge's phenomenon
%    t_node: node points
%    t_fine: fine grid
% return:
%    p_fine: the value of the interpolated polynomial at t_fine
%    p_nodes: the value of the interpolated polynomial at t_nodes

p_fine = 0; p_nodes = 0;
for kk = 1:length(t_nodes)
    p_fine = p_fine+coeff(kk,t_fine)*f(t_nodes(kk));
    p_nodes = p_nodes+coeff(kk,t_nodes)*f(t_nodes(kk));
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
end
