%Author:    Dimme de Groot
%Date:      April 2024
%Descr:     This function implemenmts the microphone spotformer
%           Inputs: 
%               audioRec        - The audio received at the microphones
%               settings        - The settings object
%               room            - The room object
%               w_mic           - The weights of the spotformer; can be left empty (i.e. w_mic = [])
%               flag_verbose    - set to true to print some info (very minimal)
%           Outputs:
%               s_out:          - The output vector containing the beamformed signal
%               w_mic:          - The computed microphone weights

function [s_out, varargout] = spotformer_microphone(audioRec, settings, room, w_mic, flag_verbose)
    if nargin == 3
        w_mic = [];
        flag_verbose = true;
    end
    if nargin == 4
        flag_verbose = true;
    end
    if isempty(w_mic)
        t = tic;
        [R_micINT, R_micTAR] = fnc_comp_CovarianceMicrophoneSpotformer(settings, room); %for microphone spotformer     
        tElap1 = toc(t);
        t = tic;
        w_mic = fnc_spotformerMic_weights(R_micTAR, R_micINT, settings)';               %Note: (transpose) conjugate!
        tElap2 = toc(t);
        if flag_verbose
            disp("Computed microphone spotformer weights in: " + num2str(tElap1+tElap2) + " seconds.")
            disp("Of this, computing the covariance matrices took: " + num2str(tElap1) + " seconds. Pretty slow! This code can likely be improved a lot")
            disp("Of this, computing the weights from the matrices took: " + num2str(tElap2) + " seconds. Pretty Quick! This code uses a generalised eigenvlaue decomposition.")
            disp(" ")
        end
    end
    flag_full_axis = settings.flag_full_axis;

    %Frame lenght, hop length, window (analysis + synthesis) 
    L1 = settings.N_t;
    R1 = L1/2;
    w_analysis = sqrthann(L1);
    w_synthesis = sqrthann(L1);
    
    %Compute beamformed playback signals
    l = 0;
    stop_flag = 0;
    s = zeros(size(audioRec,1)+settings.N_pad, room.Nr);    
    t = tic;
    while stop_flag == 0
        try
            Y_block = fft(w_analysis.*audioRec(l*R1+1:l*R1+L1,:), settings.N_t+settings.N_pad);  %FFT of single frame, zeropadded
            if flag_full_axis           
                S_block = w_mic.*Y_block;                                               %multiply with spotformerweights
                s_block = ifft(S_block);
            else
                S_block = w_mic.*Y_block(1:L1+1,:);                                     %idem
                s_block = ifft([S_block; zeros(settings.N_pad-1,room.Nr)], 'symmetric');
            end
            s(l*R1+1:l*R1+L1, :) = s(l*R1+1:l*R1+L1, :) + w_synthesis.*s_block(1:L1,:); %append frame to full signal using the windows.
        catch ME
            if strncmp("Index", ME.message, 5)
                if flag_verbose
                    disp("Crashed at frame " + num2str(l) + " v.d. approx " + num2str(floor(size(audioRec,1)/R1 - 1)))
                    disp(" ")
                end
                stop_flag = 1;
            else    
                rethrow(ME)
            end
        end
        l = l+1;
    end
    %combine the spotformed microphone signals into a single microphone signal
    s_out = sum(s,2);   
    s_out = s_out(1:size(audioRec,1));
    if flag_verbose
        tElap = toc(t);
        disp("Computing the output at the mic spotformer took an additional " + num2str(tElap) + " seconds!")
        disp(" ")
    end
    varargout{1} = w_mic;
end
