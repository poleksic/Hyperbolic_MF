function [ indices ] = OffTargetTestIndicesSet( MTX, nfolds, labind, numl )
    indices = cell(1, nfolds);

    [nrows, ncols] = size(MTX);
    nelem = nrows * ncols;
    
    if mod(nelem,nfolds) == 0
        binsize = nelem/nfolds - 1;
    else
        binsize = floor(nelem/nfolds);
    end
    
    s = RandStream('mt19937ar','Seed', sum(100*(1+labind/numl)*clock));
    RandStream.setGlobalStream(s);
    
    Perm = randperm(nelem);
    for t = 1:nfolds
        indices{t} = sort(Perm((t-1)*binsize+1:t*binsize));
        if t <= nelem-binsize*nfolds
            indices{t} = sort([indices{t} Perm(t+binsize*nfolds)]);
        end
    end
end

