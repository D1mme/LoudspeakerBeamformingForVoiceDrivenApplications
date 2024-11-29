%Author:    Dimme de Groot
%Date:      March 2024
%Descr:     Implements square root hanning window. See, for example, [1]
%       
%[1] Smith, J.O. Spectral Audio Signal Processing, http://ccrma.stanford.edu/~jos/sasp/, online book, 2011 edition.

function h = sqrthann(Lw)
    n = (0:1:Lw-1).';
    m = n-Lw/2;
    h = sqrt(ones(Lw, 1).*(0.5+0.5*cos(2*pi/Lw * m)));
end