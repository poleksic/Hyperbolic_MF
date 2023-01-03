function K = dnilmf_fastKF(sim1, sim2, nNeig, nIter)
    % same as in dnilmf

    [m,~] = size(sim1);

	S1 = sim1;
	S2 = sim2;

    S1 = S1 ./ sum(S1,2);
    S1 = (S1 + S1') / 2;
    
    S2 = S2 ./ sum(S2,2);
    S2 = (S2 + S2') / 2;
    
    S1_0 = S1;
    S2_0 = S2;
    
    [S1_sort, ~] = sort(S1,2,'descend');
    [S2_sort, ~] = sort(S2,2,'descend');
    
    for i = 1:m
        temp = S1(i,:);
        temp(temp < S1_sort(i,nNeig)) = 0;
        S1(i,:) = temp;
        
        temp = S2(i,:);
        temp(temp < S2_sort(i,nNeig)) = 0;
        S2(i,:) = temp;
    end
    
    S1 = S1 ./ sum(S1,2);    
    S2 = S2 ./ sum(S2,2);
    
    % Fusion
    for ii=1:nIter
       S1_next = S1 * S2_0 * S1';
       S2_next = S2 * S1_0 * S2';
       
       S1_0 = S1_next + eye(m);
       S2_0 = S2_next + eye(m);
       
       S1_0 = (S1_0 + S1_0')/2;
       S2_0 = (S2_0 + S2_0')/2;
    end
    
    K = (S1_0 + S2_0) / 2;
	K = K ./ sum(K, 2);
	K = (K + K' + eye(m)) / 2;    
end