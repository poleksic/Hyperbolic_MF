% R = [1,0,1;
%     0,0,0;
%     1,0,0]
% 
% N = [1, 0.2, 0.0;
%     0.2, 1, 0.3;
%     0.0, 0.3, 1.0]
% 
% NR = KNN_right(R,N,2)

% R = [1,0,1;
%     0,1,1]
% 
% N = [1, 0.2, 0.0;
%     0.2, 1, 0.3;]
% 
% NR = h_product(R,N)


m = 3;
small_positive = 0.000000000001;
U = [1,1,sqrt(3);
    2,2,sqrt(9);
    3,3,sqrt(19)]

S = [1.0,0.2,0.4;
    0.2,1.0,0.3;
    0.4,0.3,1.0]

    A = zeros(m,m);
    HP = h_product(U,U);
    for i=1:m
        A(i,i) = S(i,i);
        for j = i+1:m
            hp = HP(i,j);
            if abs(hp + 1) < small_positive
                A(i,j) = S(i,j);
            else
                A(i,j) = S(i,j)*acosh(-1*hp) / sqrt(hp^2 - 1);  
            end
            A(j,i) = A(i,j);
        end
    end
A

B = S.*acosh(-1*HP)./ sqrt(HP.*HP -ones(m,m));
B(1:m+1:end) = diag(S);
B