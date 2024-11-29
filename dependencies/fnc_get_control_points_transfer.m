%Author:    Dimme de Groot
%Date:      April 2024
%Descr:     takes as input the room, settings and the number of control points/per dim (Ntotal = Ncontrol^3).
%           returns the transfer to each of the control points and optionally also the control point locations

function [h_train, varargout] = fnc_get_control_points_transfer(room, settings)
    coorTrain = fnc_get_control_points(room.P, settings.Ncontrol);

    Ntrain = size(coorTrain, 2);
    h_train = zeros(settings.N_t, room.Ns, Ntrain);

    mtype = 'omnidirectional';  % Type of microphone
    order = -1;                 % Reflection order
    dim = 3;                    % Room dimension
    orientation = 0;            % Microphone orientation (rad)
    hp_filter = 0;              % disable high-pass filter

    %transfer to training points: direct path
    for i=1:room.Ns
        for j=1:Ntrain
            evalc("h_train(:,i,j) = rir_generator(settings.c, settings.fs, room.S(i,:), coorTrain(:,j)',room.L,[0 0 0 0 0 0], settings.N_t, mtype, order, dim, orientation, hp_filter);");
            dist(i,j) = norm(room.S(i,:)-coorTrain(:,j)');
        end
    end
    varargout{1} = coorTrain;
    varargout{2} = dist;
end