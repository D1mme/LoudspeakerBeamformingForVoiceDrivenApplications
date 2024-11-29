%Author:    Dimme de Groot
%Date:      July 2024
%Descr:     Compute the STOI value of an audio signal

function STOI = fnc_comp_STOI(clean, degraded, fs)
    [clean_a, degraded_a] = alignsignals(clean, degraded, Truncate=true);
    if length(clean_a) > length(degraded_a)
        clean_a = clean_a(1:length(degraded_a));
    else 
        degraded_a = degraded_a(1:length(clean_a));
    end
    STOI = stoi(degraded_a, clean_a, fs);
end