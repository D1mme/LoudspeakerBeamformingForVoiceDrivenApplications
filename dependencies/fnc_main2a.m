%Author:    Dimme de Groot      
%Date:      April 2024
%Descr:     Compute necessary parts for loudspeaker spotforming
%           [Wzp, L_loud, h_hat_train_diag, varargout] = fnc_main2a(settings, room, flag_full_axis)  
%                              
%Zero padded DFT matrix, decomposistion R_loud = Lloud*Lloud', the transfer to the training points h_hat_train_diag. 

function [Wzp, L_loud, h_hat_train_diag, h_hat_train, varargout] = fnc_main2a(settings, room, flag_full_axis)       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute required covariance matrices %                                                                  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    R_loud = fnc_comp_CovarianceLoudspeakerSpotformer(settings, room, flag_full_axis);  %for loudspeaker spotformer                                                         
    [L_loud, ~] = fnc_preprocess_R(R_loud, 0, false); %Decompose R_loud ~= L_loud'*Lloud
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute and store relevant stuff for loudspeaker spotformer %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [h_train, coorTrain, dist] = fnc_get_control_points_transfer(room, settings);       %Compute transfers to control points
    h_train = h_train*(4*pi*min(dist, [], 'all'));                                      %normalise based on shortest distance to control point
    h_hat_train = fft(h_train, settings.N_t+settings.N_pad);                            %frequency domain version
    
    %Now write it in a manner which is handy during loudspeaker spotforming
    N_k = size(L_loud,2);                                                   %number of freuqency bins
    h_hat_train_diag = zeros(N_k, room.Ns*N_k, settings.Ncontrol);          %we can put h_hat_train as a diagonal matrix for ease later
    for nc = 1:settings.Ncontrol
        for ns = 1:room.Ns
            h_hat_train_diag(:,(ns-1)*N_k+1:ns*N_k,nc) = sparse(diag(h_hat_train(1:N_k,ns,nc))); 
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Zero padded (and truncated) DFT matrix %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    W = dftmtx(settings.N_t+settings.N_pad);
    ZP = [eye(settings.N_t); zeros(settings.N_pad, settings.N_t)];
    Wzp = W*ZP;
    if ~flag_full_axis
        Wzp = [eye(N_k),zeros(N_k, settings.N_t+settings.N_pad-N_k)]*Wzp;   %only keep first N_k samples of zero-padded dft
    end

    % Variable outputs
    varargout{1} = h_train;
    varargout{2} = R_loud;
    varargout{3} = coorTrain;
end