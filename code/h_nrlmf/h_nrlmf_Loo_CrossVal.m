function [AUC_AVG, AUPR_AVG, PREC10_AVG, PREC100_AVG, AP_AVG] = h_nrlmf_Loo_CrossVal(R,M,N,J,c,theta,pseudo_std,alphaU,alphaV,lM,lN,eta,iter,rnk,cutoff)  

    [m n] = size(R);
        
    [DM nM]= GetDiag(M,J);
    [DN nN]= GetDiag(N,J);
    DMM = DM-nM;
    DNN = DN-nN;        
    
    AUC_AVG = 0;
    AUPR_AVG = 0;
    PREC10_AVG = 0;
    PREC100_AVG = 0;
    AP_AVG = 0;
    
    AUC = zeros(m,1);
    AUPR = zeros(m,1);
    PREC10 = zeros(m,1);
    PREC100 = zeros(m,1);
    AP = zeros(m,1);
    count_me = zeros(m,1);

    indices = ColdStartTestIndicesSet( R, m, 'CSROWS', 1, 1);

    for t=1:m
        GOLD = R(indices{t}); GOLD = GOLD(:);
        if(sum(GOLD) == 0)
            fprintf('failed row=%d\n', t);
            continue;
        else
            count_me(t) = 1;
            TEST_MTX = R;
            TEST_MTX(indices{t}) = 0;

            [~,ExcludedCols] = find(sum(TEST_MTX,1)==0);
            [ExcludedRows,~] = find(sum(TEST_MTX,2)==0);

            [SU, SV] = h_nrlmf_WeightImputeLogFactorization(TEST_MTX,DMM,DNN,c,theta,pseudo_std,alphaU,alphaV,lM,lN,eta,iter,rnk,cutoff);

            SU = h_WeightedProfile(SU, M, ExcludedRows, J);
            SV = h_WeightedProfile(SV, N, ExcludedCols, J);

            EXC = h_Get_Nrlmf_Prob(SU,SV,theta);
            EXC = EXC(indices{t});
            EXC = EXC(:);

            AUC(t) = ComputeAccuracyMeasures(GOLD,EXC,'AUC',10);
            AUPR(t) = ComputeAccuracyMeasures(GOLD,EXC,'AUPR',10);
            PREC10(t) = ComputeAccuracyMeasures(GOLD,EXC,'PREC',10);
            PREC100(t) = ComputeAccuracyMeasures(GOLD,EXC,'PREC',100);
            AP(t) = AvgPrec(GOLD,EXC,10);

% fprintf('ac=%.2f apr=%.2f prec10=%.2f\n',ac,apr,prec10);

        end
    end
    
    for t=1:m
        AUC_AVG = AUC_AVG + AUC(t);
        AUPR_AVG = AUPR_AVG + AUPR(t);
        PREC10_AVG = PREC10_AVG + PREC10(t);
        PREC100_AVG = PREC100_AVG + PREC100(t);
        AP_AVG = AP_AVG + AP(t);
    end
    AUC_AVG = AUC_AVG / sum(count_me);
    AUPR_AVG = AUPR_AVG / sum(count_me);
    PREC10_AVG = PREC10_AVG / sum(count_me);
    PREC100_AVG = PREC100_AVG / sum(count_me);
    AP_AVG = AP_AVG / sum(count_me);
end

