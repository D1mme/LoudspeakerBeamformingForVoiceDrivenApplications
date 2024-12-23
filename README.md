# Loudspeaker Beamforming for Voice Driven Applications
This repository contains the MATLAB code corresponding to the submitted paper "Loudspeaker beamforming to enhance speech recognition performance of voice driven applications" (Accepted for ICASSP 2025, link will be added once available.)

## Usage
Run `example.m`. This sets:
- the `room` object: where are the loudspeakers, the listener, etc.
- the `settings` object: this object controls the settings, such as the number of integration points in the quadrature method.
- the `par_meas` object: this object is used for the perceptual measure. See also the corresponding [paper](https://doi.org/10.1155/ASP.2005.1292) by Van de Par et al. and [my code](https://github.com/D1mme/Par-measure)

After setting the objects the loudspeaker playback signals are computed by calling the loudspeaker spotformer.

While not part of the loudspeaker spotformer itself, the performance of the loudspeaker spotformer is evaluated by considering:
- The audio quality in the neighbourhood of the listener, measured using ViSQOLAudio;
- The energy reduction achieved in the neighbourhood of the listener (compared to the reference playback file)
- The intelligibility at the output of the microphones. The microphone algorithm is either a single microphone (nearest neighbour), a microphone spotformer or an MPDR/MVDR beamformer.

Separate code for the beamformers can be found [here](https://github.com/D1mme/MPDR-beamformer) (MPDR/MVDR) and [here](https://github.com/D1mme/microphone_spotformer) (microphone spotformer).

## Requirements 
- [CVX](https://cvxr.com/cvx/) for MATLAB (I used Version 2.2, build 1148). 
  - I used [MOSEK](https://www.mosek.com/) (Version 9.1.9) as solver, but I expect other solvers to work as well.
- MATLABs [signal processing toolbox](https://www.mathworks.com/products/signal.html)
- MATLABs [audio toolbox](https://www.mathworks.com/products/audio.html)
- The [room-impulse response generator](https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator) by E. Habets. A version compiled for my system (Ubuntu) is provided.

The code was tested on Matlab R2024a on Ubuntu 23.10. 

## Notes
None at the moment :)

## Licensing
The examples make use of the [room-impulse response generator](https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator) from E. Habets (MIT license). You might need to compile this for your system.
The sound excerpt is taken from the movie ['Sprite Fight'](https://studio.blender.org/films/sprite-fright/) by Blender Studio (Creative Commons Attribution 1.0 License). The numerical integration is performed using the [Fast Clenshaw-Curtis quadrature](https://www.mathworks.com/matlabcentral/fileexchange/6911-fast-clenshaw-curtis-quadrature) by G. von Winckel, published under a permissive license. The voice commands are obtained from the [LibriSpeech ASR corpus](https://www.openslr.org/12/) and published under a CC-BY-SA 4.0 license. All licenses can be found in the corresponding folder.

The remainder of the code is published under an MIT license.

## Contact
If you find any bugs, have questions or have other comments, please contact ddegroot@tudelft.nl
