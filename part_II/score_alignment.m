function score = score_alignment(s1, s2, s2_offset)
    % count how many "same" terms there are
    score = sum(s1(1 + s2_offset : end) == s2(1: end - s2_offset));
    
    % turn that into a percentage of the overlap
    score = 100 * score / size(s2(1:end-s2_offset), 1);
end