function [AUC_AVG, AUC_STD, AUPR_AVG, AUPR_STD, PREC10_AVG, PREC10_STD, PREC100_AVG, PREC100_STD, AP_AVG, AP_STD] = ...
    h_ngn_CrossVal(R,M,N,gamma,J,c,theta,pseudo_std,alphaU,alphaV,lM,lN,eta,iter,rnk,folds,rounds,CV,cutoff)

    [m,n]=size(R);

    AUC_AVG = 0;
    AUPR_AVG = 0;
    PREC10_AVG = 0;
    PREC100_AVG = 0;
    AP_AVG = 0;
    
    AUC = zeros(rounds,1);
    AUPR = zeros(rounds,1);
    PREC10 = zeros(rounds,1);
    PREC100 = zeros(rounds,1);
    AP = zeros(rounds,1);

    S = (M + M') / 2;
    S = MakePositiveDefinite(S);
    T = (N + N') / 2;
    T = MakePositiveDefinite(T);
    
    for j=1:rounds
        indices = getCVindices(CV, R, folds, j, rounds);

        auc_sum = 0;
        apr_sum = 0;
        prec10_sum = 0;
        prec100_sum = 0;
        ap_sum = 0;
        success_folds = 0;
        t = 1;
        while t <= folds
            GOLD = R(indices{t}); GOLD = GOLD(:);
            if(sum(GOLD) == 0)
                fprintf('failed j=%d t=%d\n', j, t);
                t = t + 1;
                continue;
            else
                success_folds = success_folds + 1;
                TEST_MTX = R;
                TEST_MTX(indices{t}) = 0;
                
                [~,ExcludedCols] = find(sum(TEST_MTX,1)==0);
                [ExcludedRows,~] = find(sum(TEST_MTX,2)==0);
                
                Yr = inferZeros(TEST_MTX, S, J);
                Yc = inferZeros(TEST_MTX', T, J);
    
                KgipD = FastKgipMat(Yr, 1);
                KgipT = FastKgipMat(Yc, 1);
                
                sd = dnilmf_fastKF(KgipD,S,J,2);
                st = dnilmf_fastKF(KgipT,T,J,2);
                
                [K_d] = NeighborSimMtx(sd, J);
                [K_t] = NeighborSimMtx(st, J);

                Y_d = Construct_Y_dt(TEST_MTX, K_d);
                Y_t = Construct_Y_dt(TEST_MTX', K_t);

                A = Construct_A(K_d, K_t, TEST_MTX, Y_d, Y_t, gamma);
                H = Construct_H(K_d, K_t, TEST_MTX);

%                [DH,nH]= GetDiag(H,J); 
%                DHH = DH-nH;
       
                [G, G_bar] = MakeDiagonalMatrix(H);
                DHH = G + G_bar - (H + H');

                [SU, SV] = h_ngn_WeightImputeLogFactorization(A,DHH,c,theta,pseudo_std,alphaU,alphaV,lM,lN,eta,iter,rnk,cutoff);

                SU_t = SU(1:n,:);
                SU_d = SU(n+1:m+n,:);

                SV_t = SV(1:n,:);
                SV_d = SV(n+1:m+n,:);

%                 SU_d = h_WeightedProfile(SU_d, S, ExcludedRows, J);
%                 SV_d = h_WeightedProfile(SV_d, S, ExcludedRows, J);
%                 SU_t = h_WeightedProfile(SU_t, T, ExcludedCols, J);
%                 SV_t = h_WeightedProfile(SV_t, T, ExcludedCols, J);

                SU_d = e_WeightedProfile(SU_d, M, ExcludedRows, J, 0.0);
                SU_t = e_WeightedProfile(SU_t, N, ExcludedCols, J, 0.0);
                SV_d = e_WeightedProfile(SV_d, M, ExcludedRows, J, 0.0);
                SV_t = e_WeightedProfile(SV_t, N, ExcludedCols, J, 0.0);
                
                SU_d(:,rnk+1) = 0;
                SV_d(:,rnk+1) = 0;
                SU_t(:,rnk+1) = 0;
                SV_t(:,rnk+1) = 0;
                for i=1:m
                    SU_d(i,rnk+1) = sqrt(dot(SU_d(i,:),SU_d(i,:)) + 1);
                end
                for i=1:m
                    SV_d(i,rnk+1) = sqrt(dot(SV_d(i,:),SV_d(i,:)) + 1);
                end
                for i=1:n
                    SU_t(i,rnk+1) = sqrt(dot(SU_t(i,:),SU_t(i,:)) + 1);
                end
                for i=1:n
                    SV_t(i,rnk+1) = sqrt(dot(SV_t(i,:),SV_t(i,:)) + 1);
                end                
                
                EXC = h_ngn_getProb(SU_t,SU_d,SV_t,SV_d,theta);
                EXC = EXC(indices{t});
                EXC = EXC(:);

                ac = ComputeAccuracyMeasures(GOLD,EXC,'AUC',10);
                apr = ComputeAccuracyMeasures(GOLD,EXC,'AUPR',10);
                prec10 = ComputeAccuracyMeasures(GOLD,EXC,'PREC',10);
                prec100 = ComputeAccuracyMeasures(GOLD,EXC,'PREC',100);
                ap = AvgPrec(GOLD,EXC,10);

%fprintf('ac=%.2f apr=%.2f prec10=%.2f  prod=%1f sd = %1f\n',AUC(t),AUPR(t),PREC10(t),mn,sd);

                auc_sum = auc_sum + ac;
                apr_sum = apr_sum + apr;
                prec10_sum = prec10_sum + prec10;
                prec100_sum = prec100_sum + prec100;
                ap_sum = ap_sum + ap;
                t = t + 1;
            end
        end
        AUC(j) = auc_sum / success_folds;
        AUPR(j) = apr_sum / success_folds;
        PREC10(j) = prec10_sum / success_folds;
        PREC100(j) = prec100_sum / success_folds;
        AP(j) = ap_sum / success_folds;
        
% fprintf('ac=%.2f apr=%.2f prec10=%.2f\n',AUC(j),AUPR(j),PREC10(j));

        AUC_AVG = AUC_AVG + AUC(j);
        AUPR_AVG = AUPR_AVG + AUPR(j);
        PREC10_AVG = PREC10_AVG + PREC10(j);
        PREC100_AVG = PREC100_AVG + PREC100(j);
        AP_AVG = AP_AVG + AP(j);
    end
    AUC_AVG = AUC_AVG / rounds;
    AUPR_AVG = AUPR_AVG / rounds;
    PREC10_AVG = PREC10_AVG / rounds;
    PREC100_AVG = PREC100_AVG / rounds;
    AP_AVG = AP_AVG / rounds;
    AUC_STD = std(AUC);
    AUPR_STD = std(AUPR);
    PREC10_STD = std(PREC10);
    PREC100_STD = std(PREC100);
    AP_STD = std(AP);
end
