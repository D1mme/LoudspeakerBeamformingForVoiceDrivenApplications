%Author:    Dimme de Groot
%Date:      July 2024
%Descr:     This script is an example of how to use the microphone and loudspeaker spotformer. 
%           (1) The first section of the code (around line 45) contains some parameters. You can set
%               - the distortion parameter dPar
%               - the SIR at the microphones which is nearest to the listener
%               - Whether you want to perform loudspeaker spotforming or not (flag_loudspeaker_spotforming)
%               - The surround sound file and voice command file. 
%           (2) After this, the settings object, room and par measure are defined. 
%           The settings, room, and par measure objects are taken the same as in the corresponding submitted paper 
%               (Loudspeaker Beamforming to Enhance Speech Recognition Performance of Voice Driven Applications, submitted for ICASSP 2025)
%           (3) Then the loudspeaker spotformer is called. This returns:
%                   - the "degraded" (spotformed) playback files
%                   - the ViSQOL values at a set number of points around the listener.         
%                   - The achieved rediction in received power at a number of points in the neighbourhood of the voice assistant
%           (4) The output of the spotformer is used to perform 
%                   -   microphone spotforming
%                   -   MVDR beamforming with oracle voice activity detection (it uses the object MPDRbeamformer for that)
%                   -   MPDR beamforming 
%                   -   Nearest Neighbour "beamforming"
%               The STOI at the output of the beamformers is given
%               
%           I couldnt get the Python ASR system to work by calling it from Matlab directly, so that is excluded in this example. 
%
%Note:      (1) The Par-measure implementation is available under an MIT license at https://github.com/D1mme/Par-measure 
%           (2) On the same github, also the MPDR, MVDR and microphone spotformer can be found (possibly newer versions)
%           (3) The voice commands are taken from the LibriSpeech database. This dataset is published under a CC BY SA 4.0 license 
%           (4) The surround sound file is modified segments of the left channel of Sprite Fight, an Open Movie pubslihed by Blender Studios, it has a Creative Commons Attribution 1.0 license. 
%           (5) We use the RIR-generator by E. Habets to simulate the RIRs. IF THE CODE GIVES AN ERROR you might want to check if you need to compile it for your system. This is very easy and explained
%               in the github. See: https://github.com/ehabets/RIR-Generator . The RIR-generator is published under an MIT license. 
%           (6) For getting the Clenshaw-Curtis quadrature weights, we use "Fast Clenshaw-Curtis Quadrature" by Greg von Winckel. This function is available from  
%               https://nl.mathworks.com/matlabcentral/fileexchange/6911-fast-clenshaw-curtis-quadrature and published under a permissive license. 
%           Locations of functions mentioned in (4) and (5):
%               (in 4) dependencies/rir_generator.mexa64  - Compiled on Linux
%               (in 5) dependencies/integrate/clenquad.m
%
%           The loudspeaker spotformer is quite slow. On my laptop it takes about 4 seconds per 20 ms frame when using Mosek. I expect it to be slower for the CVX default supplied solvers. 

clear all
close all

addpath dependencies/integrate/
addpath dependencies/performance_measures/
addpath dependencies/
addpath rir_generator/
addpath clenquad/

%% The reference audio file, the distortion parameter, and whether to use only microphone spotforming or microphone + loudspeaker spotforming
dPar = 5;                                               % maximum allowable distortion
surround_file = "Surround Audio/sample_16kHz.wav";     % reference surround audio (5 channels since we use 5 loudspeakers.)
voice_command = "Voice Commands/vc2.flac";              % voice command
SIR = 0;                                                % [dB], SIR of the voice signal wrt the surround signal as measured at the nearest neighbour
flag_loudspeaker_spotforming = true;                   % True to perform loudspeaker spotforming, false to use reference signal directly          

%% Now define some flags, create a settings object, a room, the Par (distortion) measure: this is the majority of the example code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings object and flags %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Some flags
flag_reverb = true;          %[-], True if in simulated reverberant environment. False otherwise (anechoic or real)
flag_real = false;           %[-], True if in real environment. False otherwise 

flag_full_axis = false;      %[-], True for full frequency axis [-Fs/2, Fs/2). False for [0, Fs/2]. True is not tested
flag_nyquist = false;        %[-], True if the nyquist bin for the microphone sptoformer should be set to zero artificially
flag_BPF = true;             %[-], True to remove low frequencies (~<100 Hz) in the summation of the cost function of the loudspeaker spotformer
flag_par_threshold = true;   %[-], True to artifically set the low frequencies (~<100 Hz) in the weighting amtrix to a large value   
flag_donut = true;           %[-], the donut flag determines if a donut (torus) shape or a gaussian is used for integrating over the voice assistant (i.e. it sets p_M of the paper)

%Other inputs: Integration
Nx = 20;                %[-], number of points over which Clenshaw-Curtis quadrature is computed
Ny = 20;                %[-], idem
Nz = 20;                %[-], idem
Nr = 20;                %[-], idem
Ntheta = 20;            %[-], idem

VAwinrad_x = 0.5;       %[m], voice-assistant (VA) window radius (3 sigma)    (loudspeaker spotformer, no donut)
VAwinrad_y = 0.5;       %[m], voice-assistant (VA) window radius (3 sigma)    (loudspeaker spotformer, no donut)
VAwinrad_z = 0.095;     %[m], voice-assistant (VA) window radius (3 sigma)    (loudspeaker spotformer, donut & no donut)
VAwinrad_r = 0.095;     %[m], voice-assistant (VA) window radius (3 sigma)    (loudspeaker spotformer, donut)
mu_r = 0.1;             %[m], mean value in radius voice-assistant (VA)       (loudspeaker spotformer, donut)                             
INTwinrad = 0.2;        %[m], interferer window radius (3 sigma)              (microphone spotformer)    
TARwinrad = 0.5;        %[m], target window radius (3 sigma)                  (microphone spotformer)    

%Other inputs: processing
sound_vel = 342;        %[m/s], speed of sound
fs = 16000;             %[Hz], sample frequency
t_frame = 0.016;        %[s], analysis window length
t_pad = 0.016;          %[s], padding window length
if flag_reverb
    rebratio = 0.004;    %term of isotropic noise
    NUMsigma2 = 10^-8;   %numerical inacurracies
    Nsigma2 = 0;         %microphone self noise       
elseif flag_real 
    rebratio = 0.004;
    NUMsigma2 = 10^-8;
    Nsigma2 = 1e-5; 
else %anechoic
    rebratio = 0;
    NUMsigma2 = 10^-8;
    Nsigma2 = 0; 
end
N_control = 9;          %[-], the number of control points around the listener

%Create settings object
settings = Settings(sound_vel,  Nx, Ny, Nz, Nr, Ntheta, fs, t_frame, t_pad, VAwinrad_x, VAwinrad_y, VAwinrad_z,...
                    VAwinrad_r, mu_r, INTwinrad, TARwinrad, rebratio, NUMsigma2, Nsigma2, N_control, ...
                    flag_reverb, flag_real, flag_full_axis, flag_nyquist, flag_BPF, flag_par_threshold, flag_donut);

%%%%%%%%%%%%%%%
% Par-measure %
%%%%%%%%%%%%%%%
Tframe = (settings.N_t+settings.N_pad)/settings.fs;     %[s], the time of the zero-padded segmentation window. Will lead to calculating the number of samples in a single window, Nframe
x_ref = 1; x_dB_ref = 70;                               %[-], [dB SPL]; the reference value in digital and physical domain
F_cal = 1000;                                           %[Hz], the calibration frequency. Note that in the report this corresponds to f_m
Ng = 64;                                                %[-], the number of gammatone filters used
par_meas = Par_measure(settings.fs, Tframe, x_ref, x_dB_ref, F_cal, Ng);       

%%%%%%%%%%%%%%%
% Create room %
%%%%%%%%%%%%%%%
%Room walls and reflection coeffs
walls = [8.1, 6.6, 3.07];
if flag_reverb
    reflec = [0.2,-0.4, -0.6, -0.6, 0.8, 0.6];          %About 0.2 s reverberation time 
elseif flag_real
    reflec = [0.2,-0.4, -0.6, -0.6, 0.8, 0.6];          %Note that in this case its not actually "real". 
else
    reflec = [0, 0, 0, 0, 0, 0];                        %No reflections: anechoic
end

%Microphone positions
Nmic = 8;
theta = linspace(0, 2*pi, Nmic+1);      %[rad], simulate circular mic array 
theta = theta(1:Nmic); theta = -theta(:);           
r = 0.1;                                %[m], radius circular array
loc_mic = [r*cos(theta), r*sin(theta), zeros(Nmic,1)]; 
loc_mic = loc_mic + [2.01, 3.33, 0.99]; %[m], microphone positions
loc_loud = [6.010, 2.019, 1.175;        %[m], FL loudspeaker position
            0.765, 2.745, 1.175;        %[m], FR loudspeaker position
            3.370, 1.630, 1.080;        %[m], FC loudspeaker location     
            5.926, 5.330, 1.080;        %[m], BL loudspeaker location 
            2.382, 5.330, 1.275;];      %[m], BR loudspeaker location 
loc_pers = [4.01, 3.905, 1.100];        %[m], person position (target)

room = Room(loc_mic, loc_loud, loc_pers, walls, reflec, settings.c);
room.plotRoom(false);
disp("Estimated T60 time room (simulated): " + num2str(room.T60) + " [s]");

%% Compute the loudspeaker spotformer (this takes a long time!)
s_ref = settings.readAudio(surround_file); %settings.readAudio checks if sample rate matches settings.fs
flag_verbose = 1;
if flag_loudspeaker_spotforming 
    [s_out, ~] = spotformer_loudspeaker(s_ref, dPar, settings, room, par_meas, flag_verbose); %perform spotforming
else
    disp(' ')
    disp('Since flag_loudspeaker_spotforming is false, we are taking the reference playback signal as output playback signal.')
    s_out = s_ref;
end

%% Compute the ViSQOL value at a number of validation points. Also compute the reduction in received energy near the microphones
disp(" ")
if flag_loudspeaker_spotforming
    x_val = randn(4,3)*0.2 + room.P;                               %Define the validation points. You can change the number of points by setting randn(Npoints,3)
    ViSQOL = fnc_comp_ViSQOL(s_ref, s_out, x_val, room, settings);  %This function can also output the reference received and degraded received signal [ViSQOL, rec_ref, rec_deg]
    x_ene = randn(20,3)*0.1 + room.Rbar;                           %Define the points at which the energy reduction is measured.  I randomly sample points in the neighbourhood of the array center
    energyRed = fnc_comp_EnergyReduction(s_ref, s_out, x_ene, room, settings);  %[dB], this function can also output the reference received and degraded received signal [energyRed, rec_ref, rec_deg]
    disp("The mean energy reduction is " + num2str(-1*mean(energyRed)) + " dB. The standard devation is " + num2str(std(energyRed, 0)) + " dB.");
end

%% Compute the signals received at the microphones. Optionally add some inaccuracies in the locations
% Variables (see line 192/193)
% Received at points neighbourghood listener
% - audioPersonRef:     received at neighhbourhood listener when reference playback signal is playing   
% - audioPersonSpot:    received at neighbourhood listener when spotformed playback signal is playing
%
% Received at microphones VA
% - AudioRecRefVoice:   voice signal + reference playback signal
% - AudioRecSpotVoice:  voice signal + loudspeaker spotformed playback signal (if flag_spotforming true)
% - AudioRecVoice:      voice signal
% - AudioRecRef:        reference playback signal
% - AudioRecSpot:       loudspeaker sptoformed playback signal
% - audio_person:       the voice command

%In case we want to add some misplacements (otherwise place appropriate comments): 
rng("shuffle")                                                           
room_act = room;                                                            %Create the "actual room" from the modelled room. Add small inaccuracies below.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMMENT THE FOLLOWING FOUR LINES IF YOU DONT WANT INACCURACIES %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rot_xy = 5/180*pi*randn(1,1);         % Rotate microphone array (1)
room_act.R = (room.R-room.Rbar)*[cos(rot_xy) -sin(rot_xy) 0; sin(rot_xy) cos(rot_xy) 0; 0 0 1]+room.Rbar; % Rotate microphone array (2)
room_act.P = room.P + randn(size(room_act.P))*0.05;                         % Misplace listener with a standard deviation of 5 cm
room_act.S = room_act.S + randn(size(room_act.S))*0.05;                     % Misplace loudspeakers with a standard deviation of 5 cm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

audio_person = settings.readAudio(voice_command); % Read the voice command 
[audioRecRefVoice, audioRecSpotVoice, audioRecVoice, audioRecRef, audioRecSpot, audio_person] = fnc_compute_received_signal_microphones(s_out, s_ref, audio_person, room_act, settings, SIR); 

%% Compute the beamforming outputs
% Variables (e.g. see line 219 to 224): 
%   {Spot, MPDR, MVDR, NN}RefVoice - 
%       {SPOT, MPDR, MVDR, NN} is the type of microphone algorithm. 
%       Ref means reference (Ref) playback instead of loudspeaker spotformed (Spot) playback, 
%       voice means including voice command 
%   {Spot, MPDR, MVDR, NN}SpotVoice - 
%       Spot means loudspeaker spotformed (Spot) playback, instead of reference (Ref) playback, 
%       voice means including voice command 
%   {Spot, MPDR, MVDR, NN}Spot -
%       Spot means loudspeaker spotformed playback. The absence of `Voice` means that no voice-command is playinh
%   {Spot, MPDR, MVDR, NN}Ref -
%       Ref means reference playback. The absence of `Voice` means that no voice-command is included
%   {Spot, MPDR, MVDR, NN}Voice -
%       Voice means that voice command is included. The absence of {Ref, Spot} means that the loudspeakers are silent

disp(' ')
disp('Computing beamformed outputs... This might take around 2 minutes since I did not code them very efficiently.')

lambda = 0.97;
regfactor = 1e-3;
MPDR = MPDRbeamformer(settings.c, settings.fs, settings.window_length_act, settings.pad_length_act, lambda, regfactor, settings.flag_full_axis, "sqrthann", "sqrthann");

%Output after microphone spotforming 
[SpotRefVoice, w_mic] = spotformer_microphone(audioRecRefVoice, settings, room);       % (1) reference playback, voice command -> also get weights  
disp('Microphone spotformer weights computed!')
SpotSpotVoice = spotformer_microphone(audioRecSpotVoice, settings, room, w_mic, 0);    % (2) spotformed playback, voice command; verbose_flag = 0, weights are precomputed
SpotRef = spotformer_microphone(audioRecRef, settings, room, w_mic, 0);                % (3) reference playback; verbose_flag = 0, weights are precomputed
SpotSpot =  spotformer_microphone(audioRecSpot, settings, room, w_mic, 0);             % (4) spotformed playback; verbose_flag = 0, weights are precomputed
SpotVoice = spotformer_microphone(audioRecVoice, settings, room, w_mic, 0);            % (5) voice command; verbose_flag = 0, weights are precomputed

disp('Mic spotformer done...')

%Output after MVDR beamforming. The clean voice data (audioRecVoice) is used for the VAD. 
MVDRRefVoice = spotformer_microphone_MPDR(MPDR, room, audioRecRefVoice, audioRecVoice);
MVDRSpotVoice = spotformer_microphone_MPDR(MPDR, room, audioRecSpotVoice, audioRecVoice);
MVDRRef = spotformer_microphone_MPDR(MPDR, room, audioRecRef, audioRecVoice);          
MVDRSpot = spotformer_microphone_MPDR(MPDR, room, audioRecSpot, audioRecVoice);
MVDRVoice = spotformer_microphone_MPDR(MPDR, room, audioRecVoice, audioRecVoice);

disp('MVDR done...')

%Output after MPDR beamforming.
MPDRRefVoice = spotformer_microphone_MPDR(MPDR, room, audioRecRefVoice);
MPDRSpotVoice = spotformer_microphone_MPDR(MPDR, room, audioRecSpotVoice);
MPDRRef = spotformer_microphone_MPDR(MPDR, room, audioRecRef);          
MPDRSpot = spotformer_microphone_MPDR(MPDR, room, audioRecSpot);
MPDRVoice = spotformer_microphone_MPDR(MPDR, room, audioRecVoice);

disp('MPDR done ...')

%Output after nearest neighbour `beamforming`: in the paper this is referred to as NM (nearest microphone).
NNRefVoice = spotformer_nearest_neighbour(audioRecRefVoice, room);
NNSpotVoice = spotformer_nearest_neighbour(audioRecSpotVoice, room);
NNRef = spotformer_nearest_neighbour(audioRecRef, room);
NNSpot = spotformer_nearest_neighbour(audioRecSpot, room);
NNVoice = spotformer_nearest_neighbour(audioRecVoice, room);

disp('Nearest Neighbour done...')
disp('Now we are computing STOI values...')

%% Now compute the short time objective intelligiblity value 
STOI_SpotRefVoice = fnc_comp_STOI(audio_person, SpotRefVoice, settings.fs);         %microphone spotforming and reference loudspeaker playback
STOI_SpotSpotVoice =  fnc_comp_STOI(audio_person, SpotSpotVoice, settings.fs);      %microphone spotforming and spotformed loudspeaker playback
STOI_SpotVoice = fnc_comp_STOI(audio_person, SpotVoice, settings.fs);               %microphone spotforming and voice command only
STOI_MVDRRefVoice = fnc_comp_STOI(audio_person, MVDRRefVoice, settings.fs);         %MVDR beamforming and reference loudspeaker playback
STOI_MVDRSpotVoice =  fnc_comp_STOI(audio_person, MVDRSpotVoice, settings.fs);      %MVDR beamforming and spotformed loudspeaker playback
STOI_MVDRVoice = fnc_comp_STOI(audio_person, MVDRVoice, settings.fs);               %MVDR beamforming and voice command only
STOI_MPDRRefVoice = fnc_comp_STOI(audio_person, MPDRRefVoice, settings.fs);         %MPDR beamforming and reference loudspeaker playback
STOI_MPDRSpotVoice =  fnc_comp_STOI(audio_person, MPDRSpotVoice, settings.fs);      %MPDR beamforming and spotformed loudspeaker playback
STOI_MPDRVoice = fnc_comp_STOI(audio_person, MPDRVoice, settings.fs);               %MPDR beamforming and voice command only
STOI_NNRefVoice = fnc_comp_STOI(audio_person, NNRefVoice, settings.fs);             %NN "beamforming" and reference loudspeaker playback
STOI_NNSpotVoice =  fnc_comp_STOI(audio_person, NNSpotVoice, settings.fs);          %NN "beamforming" and spotformed loudspeaker playback
STOI_NNVoice = fnc_comp_STOI(audio_person, NNVoice, settings.fs);                   %NN "beamforming" and voice command only
disp("The STOI values for the microphone spotformer are (higher is better, 0 to 1): ")
disp("STOI for reference playback:  MicSpot: " + num2str(STOI_SpotRefVoice) + ". MVDR: " + num2str(STOI_MVDRRefVoice)...
        + ". MPDR: "  + num2str(STOI_MPDRRefVoice) + ". NN: " + num2str(STOI_NNRefVoice))
disp("STOI for spotformed playback (if you used loudspeaker spotforming): MicSpot: " + num2str(STOI_SpotSpotVoice) + ". MVDR: " + num2str(STOI_MVDRSpotVoice)...
        + ". MPDR: "  + num2str(STOI_MPDRSpotVoice) + ". NN: " + num2str(STOI_NNSpotVoice))

%% In case you want to hear the outputs, uncomment the following lines
disp("In case you want to hear the outputs of the different beamformers, uncomment the lines below line " +  num2str(dbstack().line+7))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNCOMMENT BELOW 50 LINES OR SO TO LISTEN TO THE OUTPUT %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SpotMax = max(abs(SpotVoice))*10; %The 10 is arbitrary but the idea is to avoid clipping. 
% MVDRMax = max(abs(MVDRVoice))*10; 
% MPDRMax = max(abs(MPDRVoice))*10; 
% NNMax = max(abs(NNVoice))*10;       
% Tpause = length(SpotRefVoice)/settings.fs;
% 
% %For the spotformer
% disp('Playing back mic spotformer outputs...')
% sound(SpotRefVoice/SpotMax, settings.fs);          
% pause(Tpause)
% sound(SpotSpotVoice/SpotMax, settings.fs);          
% pause(Tpause)
% sound(SpotRef/SpotMax, settings.fs);          
% pause(Tpause)
% sound(SpotSpot/SpotMax, settings.fs);          
% pause(Tpause)
% sound(SpotVoice/SpotMax, settings.fs);          
% pause(Tpause)
% 
% %For the MVDR
% disp('Playing back MVDR beamformer outputs...')
% sound(MVDRRefVoice/MVDRMax, settings.fs);          
% pause(Tpause)
% sound(MVDRSpotVoice/MVDRMax, settings.fs);          
% pause(Tpause)
% sound(MVDRRef/MVDRMax, settings.fs);          
% pause(Tpause)
% sound(MVDRSpot/MVDRMax, settings.fs);          
% pause(Tpause)
% sound(MVDRVoice/MVDRMax, settings.fs);          
% pause(Tpause)
% 
% %For the MPDR
% disp('Playing back MPDR beamformer outputs...')
% sound(MPDRRefVoice/MPDRMax, settings.fs);          
% pause(Tpause)
% sound(MPDRSpotVoice/MPDRMax, settings.fs);          
% pause(Tpause)
% sound(MPDRRef/MPDRMax, settings.fs);          
% pause(Tpause)
% sound(MPDRSpot/MPDRMax, settings.fs);          
% pause(Tpause)
% sound(MPDRVoice/MPDRMax, settings.fs);          
% pause(Tpause)
% 
% %For the NN
% disp('Playing back NN outputs...')
% soundsc(NNRefVoice/NNMax, settings.fs);          
% pause(Tpause)
% soundsc(NNSpotVoice/NNMax, settings.fs);          
% pause(Tpause)
% soundsc(NNRef/NNMax, settings.fs);          
% pause(Tpause)
% soundsc(NNSpot/NNMax, settings.fs);          
% pause(Tpause)
% soundsc(NNVoice/NNMax, settings.fs);         

%% 
rmpath dependencies/integrate/
rmpath dependencies/performance_measures/
rmpath dependencies/
rmpath rir_generator/
rmpath clenquad/

