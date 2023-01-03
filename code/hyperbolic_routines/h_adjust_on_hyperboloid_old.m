% Adjust the values of U(row_no,ind)

function U = h_adjust_on_hyperboloid_old(U,ind,row_no)
    [m,n] = size(U);
    U(row_no,ind) = 0;

    for i=1:length(row_no)
        if ind < n
            U(row_no(i),ind) = sqrt(2 * U(row_no(i),n)*U(row_no(i),n) - dot(U(row_no(i),:),U(row_no(i),:)) - 1);
        else
            U(row_no(i),ind) = sqrt(dot(U(row_no(i),:),U(row_no(i),:)) + 1);
        end
    end

end