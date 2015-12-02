function aligned_psi_source = construct_alignment_psi(temp, psi_source, l_shotgun, m)
    % temp: (10, 3): (num_samples, [score, other, offset])
    % psi_source: transition prob in markov chain
    
    aligned_psi_source = zeros(m,m,4,4);
    
    % Chain
    for i=1:m-1
        aligned_psi_source(i,i+1,:,:)=psi_source;
        aligned_psi_source(i+1,i,:,:)=psi_source';
    end
    % The one-off diagonals should be nonzero by now
  
    % Arcs
    num_samples = size(temp,1);
    for curr=1:num_samples
        score = temp(curr, 1);
        other = temp(curr, 2);
        offset = temp(curr, 3);
        
        len_match = l_shotgun - offset;
        
        curr_st = (curr - 1)*l_shotgun + offset;
%         curr_fin = curr_st + len_match;
        
        o_st = (other - 1)*l_shotgun+1;
%         o_fin = o_st + len_match;
        
        % Construct edge potential
        self_trans = score;
        other_trans = (1-score)/3;
        self_trans_pot = eye(4) * self_trans;
        other_trans_pot = (1-eye(4)) * other_trans;
        edge_pot = self_trans_pot + other_trans_pot;
        
        for match_idx=0:len_match-1
            aligned_psi_source(curr_st+match_idx, o_st + match_idx, :, :) = edge_pot;
            aligned_psi_source(o_st + match_idx, curr_st+match_idx, :, :) = edge_pot;
        end
    end
end