function [sc] = PrecTop(GOLDARRAY,EXCARRAY,numtop_prec)
    [B, ind] = sort(EXCARRAY(:),'descend');
    FOLLOWGOLD = GOLDARRAY(ind);
    tot_pos = sum(GOLDARRAY);
    adj_num = min(numtop_prec, tot_pos);
    sc = full(sum(FOLLOWGOLD(1:adj_num))/adj_num);
end