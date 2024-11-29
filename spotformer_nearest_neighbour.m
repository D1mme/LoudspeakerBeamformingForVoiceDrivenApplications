%Author:    Dimme de Groot
%Date:      May 2024
%Descr:     The nearest microphone algorithm: output the microphone signal which is nearest to the user 

function [s_out, varargout] = spotformer_nearest_neighbour(s_in, room)
    [L,Nmic] = size(s_in);
    if Nmic > L
        disp("In nearest_neighbour beamformer: your input signal should be L x Nmic, with L the audio length")
    end
    [~, NN] = min(vecnorm(room.R-room.P,2,2));
    s_out = s_in(:,NN);
    varargout{1}= NN;
end