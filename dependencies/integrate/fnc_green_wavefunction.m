%Author:    Dimme de Groot
%Date:      May 2024
%Descr:     This function evaluates the greens function solution to the (acoustic) wave equation in the frequency domain.
%               g(x_s, x_r, k) = exp(-1j * k * norm(x_s-x_r)) / (4*pi*norm(x_s-x_r)), with k = f/c the wavenumber, f the frequency in hertz, c the wave velocity
%           Inputs: 
%               x_s a 3x1 real vector giving the source location
%               x_r a 3xN real vector giving the receiver location
%               k a real scalar giving the wavenumber
%           Outputs:
%               out a 1xN complex vector. The i'th element contains the greens function evaluated for the i'th coordiante of x_r (i.e. x_r(:,i));
%               
%           Note that the dimensions of the inputs are important!

function out = fnc_green_wavefunction(x_s, x_r, k)    
    dist = vecnorm(x_s-x_r);
    out = 1./(4*pi*dist).*exp(-1j*k*dist);
end