function [score] = Jaccard(x,y)
    dt = dot(x,y);
    if dt == 0
        score = 0;
    else
        score = dt / sum (or(x,y));
    end
end

