%Author:    Dimme de Groot
%Date:      March 2024
%Descr:     Implements Eq. (18) of [1], with N_ss(omega) = 1 for all omega(!!!). So, N_ss(omega) needs to be incorporated later
%
%[1] Martinez et al. A robust region-based near-field beamformer. ICASSP 2015.
 
function R_iso = fnc_comp_Risotropic(k_ax, x)
    N = size(x,1);          %number of positions considered
    
    R_iso = zeros(N, N, length(k_ax));
    dist = distcalc(x,x);   %norm(x1-x2), norm(x1-x2), etc.

    for i = 1:length(k_ax)
        k = k_ax(i);
        R_iso(:,:,i) = sinc(k*dist);
    end
end