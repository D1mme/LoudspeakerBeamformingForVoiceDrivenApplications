%Author:        Dimme
%Date:          July 2024
%Description:   This function is used to compute the ViSQOL values at locations coor.
%Inputs:        - ref:  [length(ref) x room.Ns] the reference playabck (at the loudspeaker!) signal.
%               - deg:  [length(deg) x room.Ns]. With length(deg)=length(ref). the degraded playback (at the loudspeaker!) signal. 
%               - coor: [nCoor x 3]. The coordinates at which the received energy reduction should be computed
%               - room: a room object 
%               - settings: a settings object
%Outputs:       - EnergyRed:    [nCoor x 1]. dB. The energy reduction when the degraded signal is playing compared to the reference. 
%                               If settings.flag_BPF is true, the reduction is computed after high pass filtering by 100 Hz
%               - varargout{1}: [length(RIR)+length(ref)-1, nCoor], the reference signals as received at the coordinates
%               - varargout{2}: [length(RIR)+length(deg)-1, nCoor], the degraded signals as received at the coordinates
%Note:          The ViSQOL value is computed on upsampled (48 kHz, as recommended by the visqol documentation) versions of varargout{1} and varargout{2}.

function [EnergyReduction, varargout] = fnc_comp_EnergyReduction(ref, deg, coor, room, settings)
    %Number of points at which to compute ViSQOL
    [nCoor,~] = size(coor);
    EnergyReduction = nan(nCoor,1);
    
    %Compute the signal received at the point when the reference and the degraded signals are playing.
    %   (1) Compute the transfer to each of the control points. For this, we treat each control point as a person.
    room.P = coor;
    room.Np = nCoor;
    [~, ~, hLoudPerson] = room.comp_transfer(settings.fs, false);
    
    %   (2) Create placeholder data for the reference received and degraded received signal
    [nRIR, ~, ~] = size(hLoudPerson);
    nRef = size(ref,1);
    nDeg = size(deg,1);
    if nRef ~= nDeg
        disp('fnc_comp_ViSQOL expected the reference and degraded signal to be of equal length.')
    end
    rec_ref = zeros(nRIR + nRef - 1, nCoor);
    rec_deg = zeros(nRIR + nDeg - 1, nCoor);   
    
    %   (3) Compute received signals 
    for i = 1:room.Ns
        for j = 1:nCoor
            rec_ref(:,j) = rec_ref(:,j) + conv(ref(:,i), hLoudPerson(:,i,j));
            rec_deg(:,j) = rec_deg(:,j) + conv(deg(:,i), hLoudPerson(:,i,j));
        end
    end
    varargout{1} = rec_ref;
    varargout{2} = rec_deg;
    
    %Compute energy received. One catch is that it is mostly interesting over the bandwidth on which the spotformer operates. 
    %   (1) Compute HPF filter if flag_BPF is true. Note that flag_BPF is a misnomer
    f = settings.fs/2;
    if settings.flag_BPF %define 100 Hz high pass filter 
        HPFfilt = fir1(128,100/f, 'high');
    else
        HPFfilt = 1;
    end
    
    %   (2) Compute energy received at each of the coordinates
    for j=1:nCoor   
        Eref = norm(conv(rec_ref(:,j), HPFfilt));
        Edeg = norm(conv(rec_deg(:,j), HPFfilt));
        EnergyReduction(j) = 20*log10(Edeg/Eref);
    end
end