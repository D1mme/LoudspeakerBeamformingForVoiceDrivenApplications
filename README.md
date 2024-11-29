# Loudspeaker Beamforming for Voice Driven Applications
This repository contains the MATLAB code corresponding to the submitted paper "Loudspeaker beamforming to enhance speech recognition performance of voice driven applications" (submitted for ICASSP 2025)

## Usage
To be updated... For now: simply run `example.m`. In the code there is more explanation.

## Requirements 
- [CVX](https://cvxr.com/cvx/) for MATLAB installed.
  - I used [MOSEK](https://www.mosek.com/) as solver, but I expect other solvers to work as well.
- MATLABs [signal processing toolbox](https://www.mathworks.com/products/signal.html)
- MATLABs [audio toolbox](https://www.mathworks.com/products/audio.html)

## Notes
To be updated...

## Licensing
The examples make use of the [room-impulse response generator](https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator) from E. Habets (MIT license). You might need to compile this for your system.
The sound excerpt is taken from the movie ['Sprite Fight'](https://studio.blender.org/films/sprite-fright/) by Blender Studio (Creative Commons Attribution 1.0 License). The numerical integration is performed using the [Fast Clenshaw-Curtis quadrature](https://www.mathworks.com/matlabcentral/fileexchange/6911-fast-clenshaw-curtis-quadrature) by G. von Winckel, published under a permissive license (see folder).

The remainder of the code is published under an MIT license.
