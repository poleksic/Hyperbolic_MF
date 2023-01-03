function [K] = FastKF(SD, KGIP, J, iter)
% NumNei is number of elements in the neighborhood

    [m,~] = size(SD);
    PD1 = NormalizeRows(SD);
    PD1 = Symmetrize(PD1);
    PD2 = NormalizeRows(KGIP);
    PD2 = Symmetrize(PD2);

    LD1 = LocalSimilarityMatrix(PD1,J);
    LD2 = LocalSimilarityMatrix(PD2,J);
    
    PD1t = PD1;
    PD2t = PD2;
    
    for i=1:iter
        PD1tp1 = LD1 * PD2t * LD1';
        PD2tp1 = LD2 * PD1t * LD2';
        PD1t = Symmetrize(PD1tp1) + eye(m);
        PD2t = Symmetrize(PD2tp1) + eye(m); 
    end
    
    K = (PD1t + PD2t) / 2;
    K = NormalizeRows(K);
    K = (K + K' + eye(m) ) / 2;
end

