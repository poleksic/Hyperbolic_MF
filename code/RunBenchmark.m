function [] = RunBenchmark(METHOD,TARG,CV,folds,rounds)

    BASE_DIR = '/home/aleksandar/';

    addpath(strcat(BASE_DIR,'HYPERBOLIC_MF/code/euclidean/'));
    addpath(strcat(BASE_DIR,'HYPERBOLIC_MF/code/h_dnilmf/'));
    addpath(strcat(BASE_DIR,'HYPERBOLIC_MF/code/h_ngn/'));
    addpath(strcat(BASE_DIR,'HYPERBOLIC_MF/code/h_nrlmf/'));
    addpath(strcat(BASE_DIR,'HYPERBOLIC_MF/code/hyperbolic_routines/'));

    addpath(strcat(BASE_DIR,'HYPERBOLIC_MF/code/block_matrices/'));

    OUT_DIR = strcat(BASE_DIR,'HYPERBOLIC_MF/results/');

    % METHOD = 'NRLMF' or 'DNILMF' or 'NGN' or 'HYP'
    % TARG = 'Nr'; % or 'Nr', 'Gpcr', 'Ion', 'Enz', 'DrugBank', or 'Epinions', 'MovieLens'
    % CV! setting --- other choices are CV2 (CSROWS) and CV3 (CSCOLS)

    P = uploadParams(METHOD, CV, TARG);

    fn = strcat(TARG,'_M08');
    load ([fn '.dat']);
    S = full(spconvert(eval(fn)));
    fn = strcat(TARG,'_N08');
    load ([fn '.dat']);
    T = full(spconvert(eval(fn)));

    if strcmp(CV,'LOO')
        TEMP = S;
        S = T;
        T = TEMP;
    end

    [m,n] = size(P.R);

    cutoff = 500; % cutoff score for gradient trimming

    fname = strcat(strcat(strcat(strcat(METHOD,'_'),TARG),'_'),CV);
    if ~strcmp(CV,'LOO')
        fname = strcat(fname,'_CV');
    end
    fname = strcat(fname,num2str(folds));
    fid = fopen(fullfile(OUT_DIR,fname),'wt');

    if strcmp(METHOD,'HNGN')

        X=cell(length(P.rnk),length(P.J),length(P.c));
        iterations=size(X);
        pr = prod(iterations);

        numworkers = pr;

        myCluster = parcluster('local');
        myCluster.NumWorkers = numworkers;
        parpool(myCluster,numworkers);

        AUC = zeros(pr,1); AUC_STD = zeros(pr,1);
        AUPR = zeros(pr,1); AUPR_STD = zeros(pr,1);
        PREC10 = zeros(pr,1); PREC10_STD = zeros(pr,1);
        PREC100 = zeros(pr,1); PREC100_STD = zeros(pr,1);
        AP = zeros(pr,1); AP_STD = zeros(pr,1);
        ix_rnk = zeros(pr,1); ix_J = zeros(pr,1); ix_c = zeros(pr,1);

        GMAX_AUC = 0; GMAX_AUPR = 0; GMAX_PREC10 = 0; GMAX_PREC100 = 0; GMAX_AP = 0;
        GMAX_AUC_STD = 0; GMAX_AUPR_STD = 0; GMAX_PREC10_STD = 0; GMAX_PREC100_STD = 0; GMAX_AP_STD = 0;

        for gm=1:length(P.gamma)
            for re=1:length(P.eta)
                for ri=1:length(P.iter)
                    for rt=1:length(P.theta)
                        for aU=1:length(P.alphaU)
                            for lm=1:length(P.lM)
                                for ln=1:length(P.lN)
                                    for ps_std=1:length(P.pseudo_std)
                                        parfor ix=1:pr
                                            [rr, j, c]=ind2sub(iterations,ix);
                                            if strcmp(CV,'LOO')
                                                [AUC(ix), AUPR(ix), PREC10(ix), PREC100(ix), AP(ix)] = ...
                                                h_ngn_Loo_CrossVal(P.R,S,T,P.gamma{gm},P.J{j},P.c{c},P.theta{rt},P.pseudo_std{ps_std},P.alphaU{aU},P.alphaU{aU},P.lM{lm},P.lN{ln},P.eta{re},P.iter{ri},P.rnk{rr},cutoff);
                                            else
                                                [AUC(ix), AUC_STD(ix), AUPR(ix), AUPR_STD(ix), PREC10(ix), PREC10_STD(ix), PREC100(ix), PREC100_STD(ix), AP(ix), AP_STD(ix)] = ...
                                                    h_ngn_CrossVal(P.R,S,T,P.gamma{gm},P.J{j},P.c{c},P.theta{rt},P.pseudo_std{ps_std},P.alphaU{aU},P.alphaU{aU},P.lM{lm},P.lN{ln},P.eta{re},P.iter{ri},P.rnk{rr},folds,rounds,CV,cutoff);
                                            end
                                            ix_rnk(ix) = P.rnk{rr}; ix_J(ix) = P.J{j}; ix_c(ix) = P.c{c};
                                        end
                                        for ix=1:pr
                                            if AUC(ix) > GMAX_AUC
                                                GMAX_AUC = AUC(ix);
                                                GMAX_AUC_STD = AUC_STD(ix);
                                            end
                                            if AUPR(ix) > GMAX_AUPR
                                                GMAX_AUPR = AUPR(ix);
                                                GMAX_AUPR_STD = AUPR_STD(ix);
                                            end
                                            if PREC10(ix) > GMAX_PREC10
                                                GMAX_PREC10 = PREC10(ix);
                                                GMAX_PREC10_STD = PREC10_STD(ix);
                                            end
                                            if PREC100(ix) > GMAX_PREC100
                                                GMAX_PREC100 = PREC100(ix);
                                                GMAX_PREC100_STD = PREC100_STD(ix);
                                            end
                                            if AP(ix) > GMAX_AP
                                                GMAX_AP = AP(ix);
                                                GMAX_AP_STD = AP_STD(ix);
                                            end
                                            if strcmp(CV,'LOO')
                                                fprintf(fid,'AUC:%.3f [%.2f] AUPR:%.3f [%.2f] PREC10:%.3f [%.2f] PREC100:%.3f [%.2f] AP:%.3f [%.2f] eta=%.4f gamma:%.3f theta:%d iter:%d rnk:%d lM:%.2f lN:%.2f J:%d c:%d\n',AUC(ix),GMAX_AUC,AUPR(ix),GMAX_AUPR,PREC10(ix),GMAX_PREC10,PREC100(ix),GMAX_PREC100,AP(ix),GMAX_AP,P.eta{re},P.gamma{gm},P.theta{rt},P.iter{ri},ix_rnk(ix),P.lM{lm},P.lN{ln},ix_J(ix),ix_c(ix));
                                            else
                                                fprintf(fid,'AUC:%.3f(%.3f) [%.2f] AUPR:%.3f(%.3f) [%.2f] PREC10:%.3f(%.3f) [%.2f] PREC100:%.3f(%.3f) [%.2f] AP:%.3f(%.3f) [%.2f] eta=%.4f gamma:%.3f theta:%d iter:%d rnk:%d lM:%.2f lN:%.2f J:%d c:%d\n',AUC(ix),AUC_STD(ix),GMAX_AUC,AUPR(ix),AUPR_STD(ix),GMAX_AUPR,PREC10(ix),PREC10_STD(ix),GMAX_PREC10,PREC100(ix),PREC100_STD(ix),GMAX_PREC100,AP(ix),AP_STD(ix),GMAX_AP,P.eta{re},P.gamma{gm},P.theta{rt},P.iter{ri},ix_rnk(ix),P.lM{lm},P.lN{ln},ix_J(ix),ix_c(ix));                                                        
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        if strcmp(CV,'LOO')
            fprintf(fid,'HNGN: AUC:%f AUPR:%f PREC10:%f PREC100:%f AP:%f\n',GMAX_AUC,GMAX_AUPR,GMAX_PREC10,GMAX_PREC100,GMAX_AP);
        else
            fprintf(fid,'HNGN: AUC:%f(%3f) AUPR:%f(%3f) PREC10:%f(%3f) PREC100:%f(%3f) AP:%f(%3f)\n',GMAX_AUC,GMAX_AUC_STD,GMAX_AUPR,GMAX_AUPR_STD,GMAX_PREC10,GMAX_PREC10_STD,GMAX_PREC100,GMAX_PREC100_STD,GMAX_AP,GMAX_AP_STD);
        end

        poolobj = gcp('nocreate');
        delete(poolobj);
    elseif strcmp(METHOD,'HNRLMF')

        X=cell(length(P.rnk),length(P.J),length(P.c));
        iterations=size(X);
        pr = prod(iterations);

        numworkers = pr;

        myCluster = parcluster('local');
        myCluster.NumWorkers = numworkers;
        parpool(myCluster,numworkers);

        AUC = zeros(pr,1); AUC_STD = zeros(pr,1);
        AUPR = zeros(pr,1); AUPR_STD = zeros(pr,1);
        PREC10 = zeros(pr,1); PREC10_STD = zeros(pr,1);
        PREC100 = zeros(pr,1); PREC100_STD = zeros(pr,1);
        AP = zeros(pr,1); AP_STD = zeros(pr,1);
        ix_rnk = zeros(pr,1); ix_J = zeros(pr,1); ix_c = zeros(pr,1);

        GMAX_AUC = 0; GMAX_AUPR = 0; GMAX_PREC10 = 0; GMAX_PREC100 = 0; GMAX_AP = 0;
        GMAX_AUC_STD = 0; GMAX_AUPR_STD = 0; GMAX_PREC10_STD = 0; GMAX_PREC100_STD = 0; GMAX_AP_STD = 0;

        for re=1:length(P.eta)
            for ri=1:length(P.iter)
                for rt=1:length(P.theta)
                    for aU=1:length(P.alphaU)
                        for lm=1:length(P.lM)
                            for ln=1:length(P.lN)
                                for ps_std=1:length(P.pseudo_std)
                                    parfor ix=1:pr
                                        [rr, j, c]=ind2sub(iterations,ix);
                                        if strcmp(CV,'LOO')
                                            [AUC(ix), AUPR(ix), PREC10(ix), PREC100(ix), AP(ix)] = ...
                                            h_nrlmf_Loo_CrossVal(P.R,S,T,P.J{j},P.c{c},P.theta{rt},P.pseudo_std{ps_std},P.alphaU{aU},P.alphaU{aU},P.lM{lm},P.lN{ln},P.eta{re},P.iter{ri},P.rnk{rr},cutoff);                                                    
                                        else 
                                            [AUC(ix), AUC_STD(ix), AUPR(ix), AUPR_STD(ix), PREC10(ix), PREC10_STD(ix), PREC100(ix), PREC100_STD(ix), AP(ix), AP_STD(ix)] = ...
                                            h_nrlmf_CrossVal(P.R,S,T,P.J{j},P.c{c},P.theta{rt},P.pseudo_std{ps_std},P.alphaU{aU},P.alphaU{aU},P.lM{lm},P.lN{ln},P.eta{re},P.iter{ri},P.rnk{rr},folds,rounds,CV,cutoff);                                                        
                                        end
                                        ix_rnk(ix) = P.rnk{rr}; ix_J(ix) = P.J{j}; ix_c(ix) = P.c{c};
                                    end
                                    for ix=1:pr
                                        if AUC(ix) > GMAX_AUC
                                            GMAX_AUC = AUC(ix);
                                            GMAX_AUC_STD = AUC_STD(ix);
                                        end
                                        if AUPR(ix) > GMAX_AUPR
                                            GMAX_AUPR = AUPR(ix);
                                            GMAX_AUPR_STD = AUPR_STD(ix);
                                        end
                                        if PREC10(ix) > GMAX_PREC10
                                            GMAX_PREC10 = PREC10(ix);
                                            GMAX_PREC10_STD = PREC10_STD(ix);
                                        end
                                        if PREC100(ix) > GMAX_PREC100
                                            GMAX_PREC100 = PREC100(ix);
                                            GMAX_PREC100_STD = PREC100_STD(ix);
                                        end
                                        if AP(ix) > GMAX_AP
                                            GMAX_AP = AP(ix);
                                            GMAX_AP_STD = AP_STD(ix);
                                        end
                                        if strcmp(CV,'LOO')
                                            fprintf(fid,'AUC:%.3f [%.2f] AUPR:%.3f [%.2f] PREC10:%.3f [%.2f] PREC100:%.3f [%.2f] AP:%.3f [%.2f] eta=%.4f theta:%d iter:%d rnk:%d lM:%.2f lN:%.2f J:%d c:%d\n',AUC(ix),GMAX_AUC,AUPR(ix),GMAX_AUPR,PREC10(ix),GMAX_PREC10,PREC100(ix),GMAX_PREC100,AP(ix),GMAX_AP,P.eta{re},P.theta{rt},P.iter{ri},ix_rnk(ix),P.lM{lm},P.lN{ln},ix_J(ix),ix_c(ix));
                                        else
                                            fprintf(fid,'AUC:%.3f(%.3f) [%.2f] AUPR:%.3f(%.3f) [%.2f] PREC10:%.3f(%.3f) [%.2f] PREC100:%.3f(%.3f) [%.2f] AP:%.3f(%.3f) [%.2f] eta=%.4f theta:%d iter:%d rnk:%d lM:%.2f lN:%.2f J:%d c:%d\n',AUC(ix),AUC_STD(ix),GMAX_AUC,AUPR(ix),AUPR_STD(ix),GMAX_AUPR,PREC10(ix),PREC10_STD(ix),GMAX_PREC10,PREC100(ix),PREC100_STD(ix),GMAX_PREC100,AP(ix),AP_STD(ix),GMAX_AP,P.eta{re},P.theta{rt},P.iter{ri},ix_rnk(ix),P.lM{lm},P.lN{ln},ix_J(ix),ix_c(ix));                                                        
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        if strcmp(CV,'LOO')
            fprintf(fid,'HNRLMF: AUC:%f AUPR:%f PREC10:%f PREC100:%f AP:%f\n',GMAX_AUC,GMAX_AUPR,GMAX_PREC10,GMAX_PREC100,GMAX_AP);
        else
            fprintf(fid,'HNRLMF: AUC:%f(%3f) AUPR:%f(%3f) PREC10:%f(%3f) PREC100:%f(%3f) AP:%f(%3f)\n',GMAX_AUC,GMAX_AUC_STD,GMAX_AUPR,GMAX_AUPR_STD,GMAX_PREC10,GMAX_PREC10_STD,GMAX_PREC100,GMAX_PREC100_STD,GMAX_AP,GMAX_AP_STD);
        end
        poolobj = gcp('nocreate');
        delete(poolobj);
    elseif strcmp(METHOD,'HDNILMF')
        X=cell(length(P.rnk),length(P.J),length(P.c));
        iterations=size(X);
        pr = prod(iterations);

        numworkers = pr;

        myCluster = parcluster('local');
        myCluster.NumWorkers = numworkers;
        parpool(myCluster,numworkers);

        AUC = zeros(pr,1); AUC_STD = zeros(pr,1);
        AUPR = zeros(pr,1); AUPR_STD = zeros(pr,1);
        PREC10 = zeros(pr,1); PREC10_STD = zeros(pr,1);
        PREC100 = zeros(pr,1); PREC100_STD = zeros(pr,1);
        AP = zeros(pr,1); AP_STD = zeros(pr,1);
        ix_rnk = zeros(pr,1); ix_J = zeros(pr,1); ix_c = zeros(pr,1);

        GMAX_AUC = 0; GMAX_AUPR = 0; GMAX_PREC10 = 0; GMAX_PREC100 = 0; GMAX_AP = 0;
        GMAX_AUC_STD = 0; GMAX_AUPR_STD = 0; GMAX_PREC10_STD = 0; GMAX_PREC100_STD = 0; GMAX_AP_STD = 0;

        for re=1:length(P.eta)
            for ri=1:length(P.iter)
                for rt=1:length(P.theta)
                    for aU=1:length(P.alphaU)
                        for a=1:length(P.alpha)
                            for b=1:length(P.beta)
                                for g=1:length(P.gamma)
                                    for ps_std=1:length(P.pseudo_std)
                                        parfor ix=1:pr
                                            [rr, j, c]=ind2sub(iterations,ix);
                                            if strcmp(CV,'LOO')
                                                [AUC(ix), AUPR(ix), PREC10(ix), PREC100(ix), AP(ix)] = ...
                                                h_dnilmf_Loo_CrossVal(P.R,S,T,P.J{j},P.c{c},P.theta{rt},P.pseudo_std{ps_std},P.alphaU{aU},P.alphaU{aU},P.alpha{a},P.beta{b},P.gamma{g},P.iter{ri},P.rnk{rr},P.eta{re},cutoff);
                                            else
                                                [AUC(ix), AUC_STD(ix), AUPR(ix), AUPR_STD(ix), PREC10(ix), PREC10_STD(ix), PREC100(ix), PREC100_STD(ix), AP(ix), AP_STD(ix)] = ...
                                                h_dnilmf_CrossVal(P.R,S,T,P.J{j},P.c{c},P.theta{rt},P.pseudo_std{ps_std},P.alphaU{aU},P.alphaU{aU},P.alpha{a},P.beta{b},P.gamma{g},P.eta{re},P.iter{ri},P.rnk{rr},folds,rounds,CV,cutoff);                                            
                                            end
                                            ix_rnk(ix) = P.rnk{rr}; ix_J(ix) = P.J{j}; ix_c(ix) = P.c{c};
                                        end
                                        for ix=1:pr
                                            if AUC(ix) > GMAX_AUC
                                                GMAX_AUC = AUC(ix);
                                                GMAX_AUC_STD = AUC_STD(ix);
                                            end
                                            if AUPR(ix) > GMAX_AUPR
                                                GMAX_AUPR = AUPR(ix);
                                                GMAX_AUPR_STD = AUPR_STD(ix);
                                            end
                                            if PREC10(ix) > GMAX_PREC10
                                                GMAX_PREC10 = PREC10(ix);
                                                GMAX_PREC10_STD = PREC10_STD(ix);
                                            end
                                            if PREC100(ix) > GMAX_PREC100
                                                GMAX_PREC100 = PREC100(ix);
                                                GMAX_PREC100_STD = PREC100_STD(ix);
                                            end
                                            if AP(ix) > GMAX_AP
                                                GMAX_AP = AP(ix);
                                                GMAX_AP_STD = AP_STD(ix);
                                            end
                                            if strcmp(CV,'LOO')
                                                fprintf(fid,'AUC:%.3f [%.2f] AUPR:%.3f [%.2f] PREC10:%.3f [%.2f] PREC100:%.3f [%.2f] AP:%.3f [%.2f] theta=%.4f eta=%.4f iter:%d rnk:%d alpha:%.4f beta:%.4f gamma:%.4f J:%d c:%d\n',AUC(ix),GMAX_AUC,AUPR(ix),GMAX_AUPR,PREC10(ix),GMAX_PREC10,PREC100(ix),GMAX_PREC100,AP(ix),GMAX_AP,P.theta{rt},P.eta{re},P.iter{ri},ix_rnk(ix),P.alpha{a},P.beta{b},P.gamma{g},ix_J(ix),ix_c(ix));
                                            else
                                                fprintf(fid,'AUC:%.3f(%.3f) [%.2f] AUPR:%.3f(%.3f) [%.2f] PREC10:%.3f(%.3f) [%.2f] PREC100:%.3f(%.3f) [%.2f] AP:%.3f(%.3f) [%.2f] theta=%.4f eta=%.4f iter:%d rnk:%d alpha:%.4f beta:%.4f gamma:%.4f J:%d c:%d\n',AUC(ix),AUC_STD(ix),GMAX_AUC,AUPR(ix),AUPR_STD(ix),GMAX_AUPR,PREC10(ix),PREC10_STD(ix),GMAX_PREC10,PREC100(ix),PREC100_STD(ix),GMAX_PREC100,AP(ix),AP_STD(ix),GMAX_AP,P.theta{rt},P.eta{re},P.iter{ri},ix_rnk(ix),P.alpha{a},P.beta{b},P.gamma{g},ix_J(ix),ix_c(ix));
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        if strcmp(CV,'LOO')
            fprintf(fid,'HDNILMF: AUC:%f AUPR:%f PREC10:%f PREC100:%f AP:%f\n',GMAX_AUC,GMAX_AUPR,GMAX_PREC10,GMAX_PREC100,GMAX_AP);
        else
            fprintf(fid,'HDNILMF: AUC:%f(%3f) AUPR:%f(%3f) PREC10:%f(%3f) PREC100:%f(%3f) AP:%f(%3f)\n',GMAX_AUC,GMAX_AUC_STD,GMAX_AUPR,GMAX_AUPR_STD,GMAX_PREC10,GMAX_PREC10_STD,GMAX_PREC100,GMAX_PREC100_STD,GMAX_AP,GMAX_AP_STD);
        end

        poolobj = gcp('nocreate');
        delete(poolobj);
    fclose(fid);

end
