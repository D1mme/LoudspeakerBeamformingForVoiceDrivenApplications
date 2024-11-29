function [L, varargout] = fnc_preprocess_R(R, tol, flag_sub_eig)
    if nargin == 1
        tol = 0;
        flag_sub_eig = false;
    end

    Nk = size(R,3); %number of covariance matrices
    for i=1:Nk
        if flag_sub_eig %In case this flag is true, we are substracting the smallest eigenpair
            [Vsub, Dsub] = eig(R(:,:,i));
            lambda_sub = diag(Dsub);   
            [lambda_sub, sortIndxsub] = sort(lambda_sub, 'ascend');
            Vsub = Vsub(:,sortIndxsub);
            R(:,:,i) = R(:,:,i) - Vsub(:,1)*lambda_sub(1)*Vsub(:,1)'; 
        end
        [L{i}, Reig(:,i)] = fnc_CholLikeDec(R(:,:,i), tol);
    end
    varargout{1} = Reig;

    function [L, lambda] = fnc_CholLikeDec(R, tol)
        %check for symmetry
        epsilon = max(abs(R-R'));
        if epsilon < 1e-16
        else
            disp('You might want to check the symmetry, continuing anyway')
        end
        R = (R+R')/2;       %make fully symmetric
    
        %perform eigenvalue decomposition and discard eigenvalues preceeding certain threshold 
        [V,D] = eig(R);  

        %We are now gonna discard the eigenvalues below the threshold
        lambda = diag(D);   
        [lambda, sortIndx] = sort(lambda, 'ascend');
        V = V(:,sortIndx);

        if tol~=0
            INDX = find(lambda > tol*lambda(end) + eps(4)); %also discard extremely small eigenvalues        
            if lambda(max(INDX(end)-1,1))<0
                disp('There might be something wrong with the positve semidefiniteness, the one-before smallest eigenvalue is negative!')
            end
        else
            INDX = 1:length(lambda);
        end
        %the decomposition is completed as L = sqrt(lambda)*eigenvector', such that L'*L = R (apart from the smallest eigenvalue)
        L = diag(sqrt(lambda(INDX)))*V(:,INDX)';
       
      end
end



