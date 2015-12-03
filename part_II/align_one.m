function [matching_strand, matching_offset, score] = align_one(shotgun_seqs, s1_index)
% find a match from shotgun_seqs for the sequence at
% shotgun_seqs(strand_index)

    s1 = shotgun_seqs(:, s1_index);
    
    matches = [];
    for s2_index = 1:size(shotgun_seqs, 2)
        
        % don't compare the strand to itself
        if s1_index == s2_index
            continue
        end
        
        s2 = shotgun_seqs(:, s2_index);
        
        for offset = 0 : size(shotgun_seqs, 1) - 15
           score = score_alignment(s1, s2, offset);
           
           % make offsets into indices
           matches = [matches; [s2_index (offset+1) score]];
%            [s2_index (offset+1) score]
        end
    end

    [~, best_match_index] = max(matches(:, 3));
    matching_strand = matches(best_match_index, 1);
    matching_offset = matches(best_match_index, 2);
    score = matches(best_match_index, 3);
end