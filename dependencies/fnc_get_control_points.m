%Author:    Dimme de Groot
%Date:      April 2024
%Descr:     this function gets as input a matrix of coordinates of size Np x 3 and a number of points/coordinate N.
%           returns a set of control points for training and testing around the reference points. 
%

function coorTrain = fnc_get_control_points(P, N)
    disp("Note: fnc_get_control_points sets randoms seed")
    rng(12345678)       %for reproducability
    if nargin == 1
        N = 3;
    end

    Np = size(P,1);
    
    for p=1:Np
        for n = 1:N
            coorTrain(:, n, p) = 0.2/3*randn(1,3)+P(p,:);
        end
    end
end