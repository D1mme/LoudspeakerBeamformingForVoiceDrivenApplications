%Author:        Dimme
%Date:          July 2024
%Description:   This function is used to compute the ViSQOL values at locations coor.
%Inputs:        - ref:  [length(ref) x room.Ns] the reference playabck (at the loudspeaker!) signal.
%               - deg:  [length(deg) x room.Ns]. With length(deg)=length(ref). the degraded playback (at the loudspeaker!) signal. 
%               - coor: [nCoor x 3]. The coordinates at which the ViSQOL values need to be computed. 
%               - room: a room object 
%               - settings: a settings object
%Outputs:       - ViSQOL: [nCoor x 1]. The ViSQOL values at the coordinates coor
%               - varargout{1}: [length(RIR)+length(ref)-1, nCoor], the reference signals as received at the coordinates
%               - varargout{2}: [length(RIR)+length(deg)-1, nCoor], the degraded signals as received at the coordinates
%Note:          The ViSQOL value is computed on upsampled (48 kHz, as recommended by the visqol documentation) versions of varargout{1} and varargout{2}.

function [ViSQOL, varargout] = fnc_comp_ViSQOL(ref, deg, coor, room, settings)
    %Number of points at which to compute ViSQOL
    [nCoor,~] = size(coor);
    ViSQOL = nan(nCoor,1);
    
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
    
    %   (4) ViSQOL for audio (not speech!) signals recommends a sampling rate of 48 kHz. 
    rec_ref48 = resample(rec_ref, 48000, settings.fs);
    rec_deg48 = resample(rec_deg, 48000, settings.fs);

    %We now compute the ViSQOL values!
    disp('Computing ViSQOL, this might take a little while')
    for j=1:nCoor
        ViSQOL(j) = visqol(rec_deg48(:,j), rec_ref48(:,j), settings.fs);
        disp("ViSQOL " + num2str(j) + "/" + num2str(nCoor) + " done. Mean Opinion Score (MOS): " + num2str(ViSQOL(j)));
    end
end