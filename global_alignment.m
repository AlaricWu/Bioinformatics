function [score,number_of_op_solu,op_align] = global_alignment(sequence1,sequence2, score_m, score_s, score_d)
    number_of_op_solu = 1;
    op_align = 1;
    global table_m 
    global table_n
    global path_list
    J = length(sequence1) + 1;
    I = length(sequence2) + 1;
    table_n = nan*ones(I,J);
    for i = 1:I
            table_n(i,1) = 1;
        end
    for j = 1:J
            table_n(1,j) = 1;
        end
    score = D_g(I,J,sequence1,sequence2, score_m, score_s, score_d);
    path_list = [];
    find_Path(I,J,[],sequence1,sequence2, score_m, score_s, score_d);
    table_m
    optimal_score = score
    number_of_op_solu = length(path_list); 
    %number_of_op_solu = table_n(I,J); %they are equivalent, this line is more efficient.
    fprintf("The number of optimal alignment calculated by recursive search path after the table is built: %d\n", number_of_op_solu);
    fprintf("The number of optimal alignment calculated by dynamic programming when the table is building: %d\n", table_n(I,J));
    for i = 1: number_of_op_solu
        [ali1,ali2] = trace_opt_ali('g',path_list{i},sequence1,sequence2)
    end
    
end

function score = D_g(i,j,sequence1,sequence2, score_m, score_s, score_d)
    global table_m 
    global path_list
    global table_n
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
    diag = table_m(i-1,j-1)+step_score;
    up = table_m(i-1,j)+ score_d;
    left = table_m(i,j-1) + score_d;
    score = max([diag,up,left]);
    table_m(i,j) = score;
    count_path_diag = 0;
    count_path_up = 0;
    count_path_left = 0;
    if(diag == score)
        count_path_diag = table_n(i-1,j-1);
    end
    if(up == score)
        count_path_up = table_n(i-1,j);
    end
    if(left == score)
        count_path_left = table_n(i,j-1);
    end
    table_n(i,j) = count_path_diag + count_path_left + count_path_up;

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

