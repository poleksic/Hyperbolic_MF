function [ap] = AvgPrec(GOLDARRAY,EXCARRAY,k)
    % average precision at k (AP@k)
    
    [B, ind] = sort(EXCARRAY(:),'descend');
    FOLLOWGOLD = GOLDARRAY(ind);
    tot_pos = sum(GOLDARRAY);
    num_relevant_at_k = min(k, tot_pos);
    
    correct = 0;
    ap = 0;
    
    for i=1:k
        % if document at rank i is relevant
        if FOLLOWGOLD(i) == 1
            correct = correct + 1;
            ap = ap + correct / i;
        end
    end
    % divide by total number of relevant documents at k
    ap = ap / num_relevant_at_k;
end