function temp = build_temp(shotgun_seqs)
    temp = zeros(10,3);
    for s1_index = 1:size(shotgun_seqs, 1)
       [matching_strand, matching_offset, score] = align_one(shotgun_seqs, s1_index);
       temp(s1_index, :) = [score matching_strand matching_offset];
    end
end