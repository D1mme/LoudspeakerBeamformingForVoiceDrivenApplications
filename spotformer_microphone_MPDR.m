function audioOut = spotformer_microphone_MPDR(MPDR, room, audioIn, audioCleanIn, threshold)
    if nargin == 3 
        audioCleanIn = zeros(size(audioIn));
        threshold = 1e-3;
    end
    if nargin == 4
        threshold = 1e-3;
    end

    %Reset PSD matrices
    MPDR.resetPSD;

    %Define look direction 
    look_location = room.P;          %The direction in which the beamformer "looks", this is defined with respect to....
    receiver_location = room.R;    %... the receiver location 
    MPDR.updateLookDirection(look_location, receiver_location);

    %Initialise
    stop_flag = false;
    l = 0;
    audioOut = zeros(size(audioIn,1),1);
    VAD = false;    %This one can only be used in case you have oracle knowledge of the voice activity. If VAD is true the covariance matrix is not updated
    flag_bias_correction = true;

    %Loop over frames
    while ~stop_flag
        try
            audioRecFrame = audioIn(l*MPDR.N_hop+1:l*MPDR.N_hop+MPDR.N_t,:);
            audioCleanInFrame = audioCleanIn(l*MPDR.N_hop+1:l*MPDR.N_hop+MPDR.N_t,:);
            if norm(audioCleanInFrame)>threshold
                VAD = true;
            else
                VAD = false;
            end
        catch ME
            stop_flag = true;
        end
    
        if ~stop_flag
            %Compute output frame and append to outut signal
            audioOutFrame = MPDR.computeOutputFrame(audioRecFrame, flag_bias_correction, VAD); 
            audioOut(l*MPDR.N_hop+1:l*MPDR.N_hop+MPDR.N_t) = audioOut(l*MPDR.N_hop+1:l*MPDR.N_hop+MPDR.N_t) + audioOutFrame;
        
            l = l+1;
        end
    end
end