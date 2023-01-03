function [ difficultin, mostdifficultin, allin] = DifficultOffTargetTestIndicesSet( MTX )
%    BASE_DIR = '/home/aleksandar'; 
%    addpath(strcat(BASE_DIR,'/HYPERBOLIC/data/yam08/'));
%    TARG = 'Ion';
%    fn = strcat(TARG,'_R08');
%    load ([fn '.dat']);
%    MTX = full(spconvert(eval(fn)));

%     MTX = [0 0 1 1;
%           0 1 0 0;
%           1 0 0 0;
%           0 0 0 0;
%           0 0 0 0]

    rng('default');
    rng(1);

    [m,n] = size(MTX);

    fracZeros = length(find(MTX == 0)) / length(find(MTX > 0));
    
    difficultin = cell(1,1);
    mostdifficultin = cell(1,1);
    allin = cell(1,1);
    
    sumRow = sum(MTX,2);
    sumCol = sum(MTX,1);

    onediffindices = [];
    onemostdiffindices = [];
    
    for row=1:length(sumRow)
        if sumRow(row) == 1
            oneind = find(MTX(row,:) > 0);
            in = sub2ind([m,n],row,oneind);
            if sumCol(oneind) > 1
                onediffindices = [onediffindices in];
            else
                onemostdiffindices = [onemostdiffindices in];
            end
        end
    end

    for col=1:length(sumCol)
        if sumCol(col) == 1
            oneind = find(MTX(:,col) > 0);
            in = sub2ind([m,n],oneind,col);
            if any(onediffindices == in)
               onemostdiffindices = [onemostdiffindices, in];
            else
                if sumRow(oneind) > 1
                    onediffindices = [onediffindices, in];
                end
            end
        end
    end

    oneallin = [onediffindices onemostdiffindices];
    
    zerodiffindices = find(MTX == 0)';
    perm = randperm(length(zerodiffindices));
    zerodiffindices = zerodiffindices(perm);
    zerodiffindices = zerodiffindices(1:length(onediffindices));

    difficultin{1} = sort([onediffindices zerodiffindices]);
            
    zeromostdiffindices = find(MTX == 0)';
    perm = randperm(length(zeromostdiffindices));
    zeromostdiffindices = zeromostdiffindices(perm);
    zeromostdiffindices = zeromostdiffindices(1:length(onemostdiffindices));

    mostdifficultin{1} = sort([onemostdiffindices zeromostdiffindices]);
            
    
    zeroallindices = find(MTX == 0)';
    perm = randperm(length(zeroallindices));
    zeroallindices = zeroallindices(perm);
    zeroallindices = zeroallindices(1:length(oneallin));

    allin{1} = sort([oneallin zeroallindices]);    
        
    
end
