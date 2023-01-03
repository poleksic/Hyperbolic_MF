function [nF] = e_WeightedProfile(F, M, Excluded, J, WP)

    m = size(M,1);

    nF = F;
    for i = 1:m
       if sum(any(i==Excluded))>0
            Row = M(i,:);       
            Row(Excluded) = -1 * realmax;
            summation = 0;
            nF(i,:) = WP * F(i,:);
            for j=1:J
                [Mx, Ix] = max(Row);
                summation = summation + Mx;
                nF(i,:) = nF(i,:) + Mx * F(Ix,:);
                Row(Ix) = -1 * realmax;
            end
            nF(i,:) = nF(i,:) / (WP + summation);
       end      
    end
end
