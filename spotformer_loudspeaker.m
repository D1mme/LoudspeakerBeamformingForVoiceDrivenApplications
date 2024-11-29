%Author:    Dimme de Groot
%Date:      April 2024
%Descr:     This function implemenmts the loudspeaker spotformer
%           Inputs: 
%               s_ref           - The reference playack signals stored as column vector of size <audio length> x Nloudspeaker
%               settings        - A settings object with the settings for the loudspeaker spotformer 
%               room            - A room object with the room considered     
%               par_meas        - The Par masking measure
%
%           Outputs:
%               s_out:          - The output matrix containing the playback signals
%               info:           - A cell containing info about the optimisation results

function [s_out, info] = spotformer_loudspeaker(s_ref, dPar, settings, room, par_meas, verbose)
    if nargin == 6
        verbose = 0;
    end
    flag_full_axis = settings.flag_full_axis;
    flag_BPF = settings.flag_BPF;
    flag_par_threshold = settings.flag_par_threshold;

    [Wzp, L_loud, h_hat_train_diag, h_hat_train] = fnc_main2a(settings, room, flag_full_axis);    

    % In case of band pass filter, weigh L_loud appropriately 
    if flag_BPF
        %A selector function for frequencies which fall in the range of about 70 to 6800 Hz;
        w = ones(settings.N_t+settings.N_pad,1); %give the speech bandwidth more importance in the weighted summation
        INDX = max(find(par_meas.freq_ax<100)); 
        w(1:INDX) = 0;
        w(end-INDX+2:end) = 0;

        for i=1:length(L_loud)
            L_loud{i} = L_loud{i}*w(i);
        end
    end

    %Settings for framing the signal
    L1 = settings.N_t;      %frame length 
    R1 = settings.N_t/2;    %hop length 
    w1 = sqrthann(L1);          %analysis window
    w2 = sqrthann(L1);          %synthesis window



    %Start processing
    N_k = size(L_loud,2);       %number of frequency bins
    stop_flag = 0;
    s_out = zeros(size(s_ref,1), room.Ns);
    l = 0;
    while stop_flag == 0
        try
            %Get frequency domain reference (per control point) and masking curve (per control point)
            [S_ref_rec_block, p_par, S_ref_block, ~] = fnc_get_blocks(w1, l, R1, L1, s_ref, h_hat_train, settings, par_meas, flag_full_axis, flag_par_threshold);
           

            disp("Computing frame " + num2str(l+1) + " v.d. approx " + num2str(floor(size(s_ref,1)/R1 - 1)))
            tic

            if verbose == 1
                cvx_begin quiet
            else
                cvx_begin
            end
            
            cvx_solver mosek 
            variable s_hat(N_k, room.Ns) complex         %frequency domain loudspeaker signals
            variable s(L1, room.Ns)                     %loudspeaker signals
            
            obj = 0;
            for i=1:N_k
                if norm(L_loud{i})~= 0
                    obj = obj + norm(L_loud{i}*s_hat(i,:).');  %cost function   
                end
            end  
            
            minimise (obj)            
            subject to
                for i=1:settings.Ncontrol
                    norm(diag(p_par(:,i))*(h_hat_train_diag(:,:,i)*reshape(s_hat, [N_k*room.Ns,1]) - S_ref_rec_block(:,i))) <= sqrt(dPar);
                end
                s_hat == Wzp*s;
            cvx_end
            

            %Add computed signal to total signal
            s_out(l*R1+1:l*R1+L1, :) = w2.*s + s_out(l*R1+1:l*R1+L1, :);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Some info which might be handy later %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if cvx_status ~= 'solved'
                disp("There might be something fishy going on")
                flag_status = "Something fishy";
            else
                flag_status = "All fine";
            end
            status(l+1) = flag_status;
            for i=1:N_k
                cost(l+1,i) = norm(L_loud{i}*s_hat(i,:).');
                cost_ref(l+1,i) = norm(L_loud{i}*S_ref_block(i,:).');
            end
            for i=1:settings.Ncontrol
                const(l+1,i) = norm(diag(p_par(:,i))*(h_hat_train_diag(:,:,i)*reshape(s_hat, [N_k*room.Ns,1]) - S_ref_rec_block(:,i)));
            end
    
        catch ME
            if strncmp("Index", ME.message, 5)
                info{1} = cost_ref;
                info{2} = cost;
                info{3} = status;
                info{4} = const;
      
                disp("Crashed at frame " + num2str(l) + " v.d. approx " + num2str(floor(size(s_ref,1)/R1 - 1)))
                stop_flag = 1;
            else    
                rethrow(ME)
            end
        end
        l = l+1;
    end
    Nframes = l-1;
end
