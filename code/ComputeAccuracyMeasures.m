function [sc] = ComputeAccuracyMeasures(GOLDARRAY,EXCARRAY,type,numtop_prec)
    if strcmp(type, 'AUPR')
        [ ~, ~, ~, sc ] = perfcurve(GOLDARRAY, EXCARRAY, 1, 'xCrit', 'reca', 'yCrit', 'prec');
    elseif strcmp(type, 'AUC')
        [ ~, ~, ~, sc ] = perfcurve(GOLDARRAY,EXCARRAY,1);
    elseif strcmp(type, 'PREC')
        sc = PrecTop(GOLDARRAY,EXCARRAY,numtop_prec);
    else
        msg = 'computeAPR: no such measure';
        error(msg);
        exit();
    end
end