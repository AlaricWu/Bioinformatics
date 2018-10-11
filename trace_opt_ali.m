function [align_sq1,align_sq2] = trace_opt_ali(g_or_l,path, sequence1,sequence2)
    [row, col]= size(path);
    align_sq1 = [];
    align_sq2 = [];
    prev_i = path(1,col);
    prev_j = path(2,col);

    for i = 2:col
        index_i = path(1,col - i + 1);
        index_j = path(2,col - i + 1);
        if(index_i - prev_i == 1)
            align_sq2 = [align_sq2,sequence2(index_i -1)];
        else
            align_sq2 = [align_sq2,' ']; 
        end
        if(index_j - prev_j == 1)
            align_sq1 = [align_sq1,sequence1(index_j -1)];
        else
            align_sq1 = [align_sq1,' ']; 
        end
        prev_i = index_i;
        prev_j = index_j;
    end
end