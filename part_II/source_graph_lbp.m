function [M_from_source, messages] = source_graph_lbp(M_to_source, messages,...
    phi_source, aligned_psi_source, r, w)
    % = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
    % [3] SOURCE GRAPH BP ITERATION
    
    % **************************************************************
    % ************ write your source graph BP code here ************
    % **************************************************************
    % (either directly in here or as function call)
    % input variables:
    %   M_to_source - msg passed from code to source graph
    %                 [size m x 4]
    %   messages - msgs in source graph that is being updated 
    %                 [size m x m x 4 x 2]
    %   phi_source - doped node potentials 
    %                [size m x 4] (given by us)
    %   aligned_psi_source - the representative edge potential matrix 
    %                representing probabilistic transitions in markov chain
    %                from i to j [size m x m x 4 x 4] 
    % output variable:
    %   M_from_source - product of incoming source graph msgs at every node
    %                   to be passed to code graph (all after msg update)
    %                   [size m x 4]
    %   messages - updated msgs in source graph
    %                 [size m x m x 4 x 2]
    
    % NOTE: psi_source for chain is not symmetrical! Special case for
    % backward msgs?
    
    m = size(phi_source,1);
    
    % LBP
    for j=1:m
        nj = find(aligned_psi_source(j,:,1,1) > 0);  % (1, num_nj) the 1,1 are arbitrarily chosen because the transition matrix has positive prob everywhere
        for xj=1:4
            for i=nj  % iterating by value
                messages(i,j,xj,w) = 0;
                for xi=1:4
                    ni = find(aligned_psi_source(i,:,1,1) > 0);  % (1, num_ni)
                    ni = ni(ni ~= j);
                    msgs = [];   % row vector (1, num_ni
                    for k=ni
                        msgs = [msgs messages(k,i,xi,r)];
                    end
                    p = prod(msgs);
                    result = phi_source(i, xi) * M_to_source(i,xi) * p * aligned_psi_source(i,j,xi,xj);
                    messages(i,j,xj,w) = messages(i,j,xj,w) + result;
                end
            end
        end
    end
    
    % Normalize LBP
    for j=1:m
        nj = find(aligned_psi_source(j,:,1,1) > 0);  % the 1,1 are arbitrarily chosen because the transition matrix has positive prob everywhere
        for i=nj
            messages(i,j,:,w) = messages(i,j,:,w)/sum(messages(i,j,:,w));
        end
    end
    
    % M_from_source
    for j=1:m
        nj = find(aligned_psi_source(j,:,1,1) > 0);  % the 1,1 are arbitrarily chosen because the transition matrix has positive prob everywhere
        for xj=1:4
            p = prod(messages(nj,j,xj,w));
            M_from_source(j,xj)=p;
        end
    end
    
    % Normalize M_from_source
    mfs_sum = sum(M_from_source,2);
    M_from_source = bsxfun(@rdivide, M_from_source, mfs_sum);
    
    
%     
%     
% 
%     % Run forward step
%     % forward(i, sip1, :) is message from i to i+1 for alphabet of s_{i+1}
%     forward(1,:,w) = (phi_source(1,:) .* M_to_source(1,:)) * psi_source;
% %     for i = 2:m-1
% %         forward(i,:,w) = (phi_source(i,:) .* M_to_source(i,:) .* (forward(i-1,:,r))) * psi_source;
% %     end
%     forward(2:m-1,:,w) = (phi_source(2:m-1,:) .* M_to_source(2:m-1,:) .* (forward(1:m-2,:,r))) * psi_source;
% 
%     
%     % Normalize forward
%     fsum = sum(forward(:,:,w),2);
%     forward(:,:,w) = bsxfun(@rdivide, forward(:,:,w), fsum);
%     
%     % Run backward step
%     % backward(i, sim1, :) is message from i to i-1 for alphabet of s_{i-1}    
%     backward(1,:,w) = (phi_source(1,:) .* M_to_source(1,:)) * psi_source';
% %     for i = m-1:-1:2
% %         backward(i,:,w) = (phi_source(i,:) .* M_to_source(i,:) .* backward(i+1,:,r)) * psi_source';
% %     end
%     backward(m-1:-1:2,:,w) = (phi_source(m-1:-1:2,:) .* M_to_source(m-1:-1:2,:) .* backward(m:-1:3,:,r)) * psi_source';
% 
% %     % Normalize backward
% %     bsum = sum(backward(:,:,w),2);
% %     backward(:,:,w) = bsxfun(@rdivide, backward(:,:,w), bsum);
% %     
%     % Compute M_from_source
%     M_from_source(1,:) = backward(2,:,w) .* phi_source(1,:);
%     for i = 2:m-1
%         M_from_source(i,:) = forward(i-1,:,w) .* backward(i+1, :, w) .* phi_source(i,:);
%     end
%     M_from_source(m,:) = forward(m-1,:,w) .* phi_source(m,:);
% 
%     % Normalize M_from_source
%     mfs_sum = sum(M_from_source,2);
%     M_from_source = bsxfun(@rdivide, M_from_source, mfs_sum);
 
end