%Author:    Dimme de Groot
%Date:      March 2024
%Descr:     computes the signal  received at the microphones. 
%               
%Outputs:   audioRec:   the received signal of size N x Nr, with N = length(longest input signal) + length(RIR) - 1
%Inputs:    interferer: the set of interferer (loudspeaker) signals. Given as a matrix of size L x Ns, with L the signal length
%           target:     the target signal, given as a cell of size Np
%           room:       the room object
%           settings:   the settings object
function [audioRecRefVoice, audioRecSpotVoice, audioRecVoice, audioRecRef, audioRecSpot, voice_command] = fnc_compute_received_signal_microphones(loud_out, loud_ref, voice_command, room, settings, SNR)
    %Compute length of each command
    nVoice = length(voice_command);
    nRef = size(loud_ref,1);
    nOut = size(loud_out,1);
    if nRef ~= nOut
        disp('expected reference playback signal and spotformed playback signal to have equal length')
    end

    %Introduce a random delay to the voice signal and make it equal length to to the loudspeaker playback signals
    nDelay = randi([1, nRef-nVoice-1],1,1);
    voice_command = [zeros(nDelay,1); voice_command; zeros(nRef - nVoice - nDelay,1)];
    nVoice = length(voice_command); %update length. 
    if nRef ~= nVoice
        disp('expected voice signal and reference signal to have equal length. If you see this there is a bug in fnc_compute_received_signal or your voice command is too long')
    end

    %Lets compute the received signals!
    %   (1) Compute RIRs
    [hLoud, hPerson, ~] = room.comp_transfer(settings.fs, false); %Get transfer from Loudspeaker to Microphones and from Person to Microphones
    [nRIR, ~, ~] = size(hLoud);

    audioRecVoice = zeros(nRIR+nVoice-1, room.Nr);
    audioRecRef = zeros(nRIR+nVoice-1, room.Nr);
    audioRecSpot = zeros(nRIR+nVoice-1, room.Nr);

    %   (2) Compute signals received from loudspeaker playback
    for nLoud = 1:room.Ns 
        for nRec = 1:room.Nr
            audioRecRef(:,nRec) = audioRecRef(:,nRec) + conv(loud_ref(:,nLoud), hLoud(:,nLoud, nRec));
            audioRecSpot(:,nRec) = audioRecSpot(:,nRec) + conv(loud_out(:,nLoud), hLoud(:,nLoud, nRec));
        end
    end

    %   (3) Compute signals received from voice command
    for nPerson = 1:room.Np 
        for nRec = 1:room.Nr
            audioRecVoice(:,nRec) = audioRecVoice(:,nRec) + conv(voice_command(:,nPerson), hPerson(:,nPerson,nRec));
        end
    end
    
    %We are now gonna rescale the voice command (measured at the microphone nearest to the person) the desired SNR
    I = find(abs(audioRecVoice(:,room.NN))>0.05*max(abs(audioRecVoice(:,room.NN))));
    S = norm(audioRecVoice(I(1):I(end), room.NN));
    N = norm(audioRecRef(I(1):I(end),room.NN));
    alpha = 10^(SNR/20) *N/S;

    audioRecVoice = alpha*audioRecVoice;
    audioRecRefVoice = audioRecRef + audioRecVoice;
    audioRecSpotVoice = audioRecSpot + audioRecVoice;
end