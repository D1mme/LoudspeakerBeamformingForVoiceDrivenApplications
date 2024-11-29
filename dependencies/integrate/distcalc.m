%Author:    Dimme de Groot
%Date:      Narch 2024
%Descr:     D = distcalc(X1, X2) computes the distances between X1 and X2 for each item in the list
%               X1 is an M1 x dim matrix with entires (x,y,z) or similar (dim the number of dimensions)
%               X2 is an M2 x dim "                                                                   " 
%               D is an M1 x M2 matrix with the distances between the locations in X1 and X2. The i'th row of the j'th columnt corresponds to L2-norm( {X1}_i - {X2}_j ) 

function D=distcalc(X1,X2)
    M1 = size(X1,1);
    M2 = size(X2,1);
    D = zeros(M1,M2);
    for i=1:M1
        D(i,:) = vecnorm(X1(i,:) - X2,2,2);
    end
end
