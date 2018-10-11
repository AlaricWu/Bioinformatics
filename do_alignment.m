function [score,number_of_op_solu,op_align] = do_alignment(g_or_l, sequence1,sequence2, score_m, score_s, score_d)
    J = length(sequence1) + 1;
    I = length(sequence2) + 1;
    global table_m 
    table_m = nan*ones(I,J);
    
    if(g_or_l == 'g')
        for i = 1:I
            table_m(i,1) = (i-1)*score_d;
        end
        for j = 1:J
            table_m(1,j) = (j-1)*score_d;
        end
    [score,number_of_op_solu,op_align] = global_alignment(sequence1,sequence2, score_m, score_s, score_d);
    
    elseif(g_or_l == 'l')
        for i = 1:I
            table_m(i,1) = 0;
        end
        for j = 1:J
            table_m(1,j) = 0;
        end
    [score,number_of_op_solu,op_align] = local_alignment(sequence1,sequence2, score_m, score_s, score_d);
    end
        
        
end

