% Adjust the values of U(row_no)

function U = h_adjust_on_hyperboloid_mtx(U,row_no)

    for i=1:length(row_no)
        U(row_no(i),:) = h_adjust_on_hyperboloid(U(row_no(i),:));
    end

end