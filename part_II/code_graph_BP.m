function [M_from_code, o_code] = code_graph_BP(M_to_code, o_code, x, H, phi_code)
    [k,n] = size(H);
    M_from_code = zeros(n, 2);
    Mhat_to_code = zeros(n, 1);
    for i = 1:n
        Mhat_to_code(i) = log(M_to_code(i, 1)) - log(M_to_code(i, 2));
    end
%     display(sum(abs(M_to_code(:, 1) - M_to_code(:, 2))))
    
    phi_hat = ones(n, 1);
    for i = 1:size(phi_hat, 1)
        phi_hat(i) = log(phi_code(i, 1)) - log(phi_code(i, 2));
    end
    
    % node to factor
%     display('node to factor')
    for i = 1:n
%         fprintf(['node' num2str(i) '\n']);
        node_pot = phi_hat(i);
        for a = 1:k 
            if H(a, i) == 0
                continue;
            end
            
            % use the node potential for this node
            msg = node_pot;
            
            % use information coming from the source graph
            msg = msg + Mhat_to_code(i);
            
            for b = 1:k
                if b == a || H(b, i) == 0
                    continue;
                end
                
                msg = msg + o_code.factor_to_node(b, i);
            end
            o_code.node_to_factor(i, a) = msg;
        
        end
    end
        
    
    % factor to node
%     display('factor to node')
    for a = 1:k
%         fprintf(['factor' num2str(a) '\n']);
        for i = 1:n
            if H(a, i) == 0
                continue;
            end
            
            msg = (-1).^(x(a));
                        
            for j = 1:n
                if i == j || H(a, j) == 0
                    continue;
                end
                msg = msg * tanh( o_code.node_to_factor(j, a) / 2 );
            end
            msg = 2 * atanh(msg);
            o_code.factor_to_node(a, i) = msg;
        end
    end
    
    % write out the M_from_code
    
    m0_tensor = 0.5 * (1 + tanh(0.5 * o_code.factor_to_node(:, :)));
    m1_tensor = 0.5 * (1 - tanh(0.5 * o_code.factor_to_node(:, :)));
    for i = 1:n
        M_from_code(i, 1) = phi_code(i, 1) .* prod(m0_tensor(find(H(:, i))',i));
        M_from_code(i, 2) = phi_code(i, 2) .* prod(m1_tensor(find(H(:, i))',i));
        M_from_code(i, :) = M_from_code(i, :) / sum(M_from_code(i, :));
        
%         m0 = (1 + tanh(
        
%         msg = 1;
%        
%         for a = 1:k
%             m0_minus_m1 = tanh(o_code.factor_to_node(a, i) / 2);
%             msg = msg * m0_minus_m1;
%         end
% %         m0_minus_m1 = tanh(msg / 2);
% 
%         m0 = (1 + m0_minus_m1) / 2;
%         m0 = m0 * phi_code(i, 1);
%         
%         m1 = 1 - m0;
%         m1 = m1 * phi_code(i, 2);
%         
%         M_from_code(i, :) = [m0, m1];
%         M_from_code(i, :) = M_from_code(i, :) / sum(M_from_code(i, :));
    end
end
