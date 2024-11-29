%Author:    Dimme de Groot
%Date:      March 2024
%Descr:     This function computes the weights of the microphone spotformer. 
%
%Inputs:    R_TAR in Nr x Nr x Nk x Np:     the covariance matrices of the target locations
%           R_INT in Nr x Nr x Nk x Nint:   the covariance matrices of the interferer locations
%           settings                        the setting object; we use settings.flag_full_axis and settings.nyquist
%Outputs:   w in Nr x Nk:                   the beamforming weights. 
%           lambda in Nk (I think)          the eigenvalues of the GEVD
%with Nr the number of microphones, Np the number of targets, Nint the number of interfereres and Nk the number of frequency bins


function [w, varargout] = fnc_spotformerMic_weights(R_TAR, R_INT, settings)
    full_axis_flag = settings.flag_full_axis;
    nyquist_flag = settings.flag_nyquist;    

    N_k = size(R_TAR,3);    %number of frequency bins

    for k=1:N_k
        R_INTk = R_INT(:,:,k);    %compute total interferer covariance by summing over the individual ones
        R_TARk = R_TAR(:,:,k);    %idem      
        [V,D] = eig(R_INTk, R_TARk);        %perform generalised eigenvalue decomposition. 
        d = diag(D);
        [min_d,i] = min(d);
        
        %We expect that the imaginary part of d is virtually zero. Similarly, its eigenvalue should be larger than 0 (postive definite)
        if(sum(abs(imag(d)))~=0)
            disp("Eigenvalues of the " + num2str(k) + "th bin are partially imaginary." ...
            + "Value (sum of absolute imagninary parts)=" + num2str(sum(abs(imag(d)))))
        end
        if min_d < 0
            disp("Numerical inaccuracy (I think): the smallest eigenvalue of bin "...
                + num2str(k) + " is below zero. Value: " +num2str(min_d) );
        end

        %Take eigenvector of length 1 corresponding to smallest eigenvalue.
        v = V(:,i);
        w(:,k) = v/norm(v);
        lambda(k) = min_d;
    end

    %Artifically set Nyquist bin to zero in case nyquist_flag = true
    if nyquist_flag 
        if full_axis_flag
            w(:,N_k/2+1) = 0;
        else
            w(:,end) = 0;
        end
    end
    varargout{1} = lambda;
end