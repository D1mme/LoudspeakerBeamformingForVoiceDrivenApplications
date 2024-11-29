%Author:    Dimme de Groot
%Date:      March 2024
%Descr:     This function is a wrapper for computing the covariance matrix R_L = R_h + R_iso + R_num of the loudspeaker spotformer
%
%Inputs:
%       - settings object
%       - room object
%       - full_axis_flag: this flag specifies whether or not to use the full frequency axis or to use conjugate symmetry when computing the PSD matrices
%Outputs:
%       - R_L the PSD matrix of the loudspeaker spotformer

function [R_L, varargout] = fnc_comp_CovarianceLoudspeakerSpotformer(settings, room, full_axis_flag)
    %general stuff
    N = [settings.Nx, settings.Ny, settings.Nz, settings.Nr, settings.Ntheta];
    fs = settings.fs;
    N_t = settings.N_t;
    N_pad = settings.N_pad;

    sigma2 = [settings.VAwinsigma2x, settings.VAwinsigma2y, settings.VAwinsigma2z, settings.VAwinsigma2r];
    x = room.S;
    
    %Define frequency-axis
    N_k = N_t + N_pad;                              %[-], total number of frequency bins (even)  
    if full_axis_flag
        k_ax = (-N_k/2:1:N_k/2-1)/N_k*2*pi*fs/room.c; 
        k_ax = fftshift(k_ax); 
    else
        k_ax = (0:1:N_k/2)/N_k*2*pi*fs/room.c;          %[rad/m], wave number -> only half of the axis [0, Fs/2] is relevant
    end
   
    mu = [room.Rbar, settings.mu_r];
    R_h = fnc_comp_Rregion(N, k_ax, x, mu, sigma2, settings.flag_donut);

    R_num = fnc_comp_Rnum(settings.NUMsigma2, room.Ns);                      
    R_iso = fnc_comp_Risotropic(k_ax, x);
    R_L = R_h + R_num;                      %I assume R_num to be added to R_h directly
    normalisation = max(pagenorm(R_L, "fro"), [], 'all');
    
    R_L = R_L./normalisation;        %Normalise frobenius norm per frequency bin
    R_L = R_L + settings.rebratio*pagenorm(R_L, "fro").*R_iso;    %Now rebratio can be added

    %in case the user wants the frequency axis back. 
    f_ax = k_ax*room.c/(2*pi);     %[rad/s], k_ax*room.c -> only half of the axis [0, Fs/2] is relevant
    varargout{1} = f_ax;
end