%Author:    Dimme de Groot
%Date:      March 2024
%Descr:     This function is a wrapper for computing the covariance matrix R_L = R_h + R_iso + R_num of the microphone spotformer
%
%Inputs:    Settings and Room objects
%
%Outputs:   R_L_INT: the interferer covariance
%           R_L_TAR: the target spotformer
%           varargout: radial frequency axis

function [R_L_INT, R_L_TAR, varargout] = fnc_comp_CovarianceMicrophoneSpotformer(settings, room) 
    %general stuff
    full_axis_flag = settings.flag_full_axis;
    N = [settings.Nx, settings.Ny, settings.Nz];
    
    fs = settings.fs;
    N_t = settings.N_t;
    N_pad = settings.N_pad;

    sigma2_INT = settings.INTwinsigma2;   
    sigma2_TAR = settings.TARwinsigma2;   
    x = room.R;
    
    %Define frequency-axis
    N_k = N_t + N_pad;                              %[-], total number of frequency bins (even)  
    if full_axis_flag
        k_ax = (-N_k/2:1:N_k/2-1)/N_k*2*pi*fs/room.c; 
        k_ax = fftshift(k_ax); 
    else
        k_ax = (0:1:N_k/2)/N_k*2*pi*fs/room.c;          %[rad/m], wave number -> only half of the axis [0, Fs/2] is relevant
    end
                     
    %Compute covariances
    R_L_INT = 0;
    for i=1:room.Ns  %each loudspeaker is an interferer
        mu = room.S(i,:);
        R_h = fnc_comp_Rregion(N, k_ax, x, mu, sigma2_INT);    %compute the covariance over the regions. Either using matlabs integral3 or gaussian-hermite(depending on matlab_flag)
        R_L_INT = R_L_INT + R_h; %We sum the contribution of each region to get the total region
    end
    R_n = fnc_comp_Rnum(settings.Nsigma2, room.Nr);
    R_iso = fnc_comp_Risotropic(k_ax, x);
    R_num = fnc_comp_Rnum(settings.NUMsigma2, room.Nr);     
    R_L_INT = R_L_INT + R_num + R_n + settings.rebratio*pagenorm(R_L_INT).*R_iso;

    R_L_TAR = 0;
    for i=1:room.Np
        mu = room.P(i,:);   %each person is a target; So far I only used one person so it might not work with more than one.
        R_h = fnc_comp_Rregion(N, k_ax, x, mu, sigma2_TAR, false);    %compute the covariance over the regions.                  
        R_L_TAR = R_L_TAR+R_h; %We sum the contribution of each region to get the total region
    end 
    R_n = fnc_comp_Rnum(settings.Nsigma2, room.Nr);
    R_num = fnc_comp_Rnum(settings.NUMsigma2, room.Nr);   
    R_L_TAR = R_L_TAR + R_num + R_n; %note: no R_iso added here

    %in case the user wants the frequency axis back. 
    w_ax = k_ax*room.c;     %[rad/s], k_ax*room.c -> only half of the axis [0, Fs/2] is relevant
    varargout{1} = w_ax;
end