function [indices] = getCVindices(CV, R, folds, j, rounds)
    if strcmp(CV, 'OFFTARGET')
        indices = OffTargetTestIndicesSet( R, folds, j, rounds);
    elseif strcmp(CV, 'COLDSTART') 
        indices = ColdStartTestIndicesSet( R, folds, CV, j, rounds);
    end
end
