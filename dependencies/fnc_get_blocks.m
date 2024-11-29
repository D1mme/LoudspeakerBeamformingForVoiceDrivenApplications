%Author:    Dimme de Groot
%Date:      03-04-204
%Descr:     just a simple function which is used to compute the par curves and the received signal at each of the control points.

function [s_hat_block, p_par, varargout] = fnc_get_blocks(w1, l, R1, L1, s, h_hat, settings, par_meas, flag_full_axis, flag_par_threshold)
    s_block = w1.*s(l*R1+1:l*R1+L1,:);
    frame_hat_block = fft(s_block, settings.N_t+settings.N_pad);    %Playback signal of frame
    s_hat_block = squeeze(sum(frame_hat_block.*h_hat,2));           %Received signal at control point, freq. domain
    s_block_par = ifft(s_hat_block);

    %Compute maskcurve
    for i=1:settings.Ncontrol
        if flag_par_threshold
            [~, ~, p_par(:, i)] = par_meas.comp_maskcurve(s_block_par(:,i), false, 100); %Set frequency bins below 100 high artificially
        else
            [~, ~, p_par(:, i)] = par_meas.comp_maskcurve(s_block_par(:,i));
        end
    end

    %half frequency axis
    if ~flag_full_axis
        N_k = (settings.N_t+settings.N_pad)/2+1;        %length halfed frequency axis
        
        p_par = p_par(1:N_k,:);                         %half p_par               
        p_par(2:N_k-1,:) = p_par(2:N_k-1,:)*sqrt(2);    %needed to keep total cost constant even when halfing frequency axis     
        s_hat_block = s_hat_block(1:N_k,:);
        frame_hat_block = frame_hat_block(1:N_k,:);
    end
    if nargout > 2
        varargout{1} = frame_hat_block;
        varargout{2} = s_block;
    end
end