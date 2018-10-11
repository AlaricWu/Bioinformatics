function [score,number_of_op_solu,op_align] = global_alignment(sequence1,sequence2, score_m, score_s, score_d)
    number_of_op_solu = 1;
    op_align = 1;
    global table_m 
    global path_list
    J = length(sequence1) + 1;
    I = length(sequence2) + 1;
    score = D_g(I,J,sequence1,sequence2, score_m, score_s, score_d);
    path_list = [];
    find_Path(I,J,[],sequence1,sequence2, score_m, score_s, score_d);
    table_m
    number_of_op_solu = length(path_list);
    for i = 1: number_of_op_solu
        [ali1,ali2] = trace_opt_ali('g',path_list{i},sequence1,sequence2)
    end
    score
end

function score = D_g(i,j,sequence1,sequence2, score_m, score_s, score_d)
    global table_m 
    if(isnan(table_m(i-1,j-1)))
        table_m(i-1,j-1) = D_g(i-1,j-1,sequence1,sequence2, score_m, score_s, score_d);
    end
    step_score = 0;
    if(sequence1(j-1) == sequence2(i-1))
        step_score = score_m;
    else
        step_score = score_s;
    end

    if(isnan(table_m(i,j-1)))
        table_m(i,j-1) = D_g(i,j-1,sequence1,sequence2, score_m, score_s, score_d);
    end
    
    if(isnan(table_m(i-1,j)))
        table_m(i-1,j) = D_g(i-1,j,sequence1,sequence2, score_m, score_s, score_d);
    end
    score = max([table_m(i-1,j-1)+step_score , table_m(i-1,j)+ score_d,table_m(i,j-1) + score_d]);
    table_m(i,j) = score;

end

function find_Path(i,j,prev_path,sequence1,sequence2, score_m, score_s, score_d)
    global path_list
    global table_m
    cur_path = [prev_path,[i;j]];
    if(i == 1 && j == 1)
        path_list{length(path_list) + 1} = cur_path;
        return;
    end
    step_score = 0;
    if(sequence1(j-1) == sequence2(i-1))
        step_score = score_m;
    else
        step_score = score_s;
    end

    diag = table_m(i-1,j-1) + step_score;
    up = table_m(i-1,j) + score_d;
    left = table_m(i,j-1) + score_d;
    if(diag == table_m(i,j))
        find_Path(i-1,j-1,cur_path,sequence1,sequence2, score_m, score_s, score_d);
    end
    if(up == table_m(i,j))
        find_Path(i-1,j,cur_path,sequence1,sequence2, score_m, score_s, score_d);
    end
    if(left == table_m(i,j))
        find_Path(i,j-1,cur_path,sequence1,sequence2, score_m, score_s, score_d);
    end
end

