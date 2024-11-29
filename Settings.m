%Author:    Dimme de Groot
%Date:      March 2024
%
%Functions: 
%   settings = Settings(<>)    
%   y = settings.readAudio(audiopath)
%   settings = settings.updSigma(INTwinrad, TARwinrad) 
classdef Settings
    properties
        c                   %[m/s], speed of sound
        Nx                   %[-], nbr of points in Gauss-Hermite quadrature
        Ny                   %[-], nbr of points in Gauss-Hermite quadrature
        Nz                   %[-], nbr of points in Gauss-Hermite/Clenshaw-curtis quadrature
        Nr                   %[-], nbr of points in Gauss-Hermite/Clenshaw-curtis quadrature
        Ntheta               %[-], nbr of points in Clenshaw-Curtis quadrature
        fs                  %[Hz],  sampling frequency

        window_length_act   %[s], actual window length
        N_t                 %[-], window length in samples      
        pad_length_act      %[s], actual padding length
        N_pad               %[-], padding length in samples

        VAwinsigma2x         %[m^2], sigma^2 window Voice Assistant  (loudspeaker spotformer)   
        VAwinsigma2y         %[m^2], sigma^2 window Voice Assistant  (loudspeaker spotformer) 
        VAwinsigma2z         %[m^2], sigma^2 window Voice Assistant  (loudspeaker spotformer) 
        VAwinsigma2r         %[m^2], sigma^2 window Voice Assistant  (loudspeaker spotformer) 
        mu_r                %[m], center of donut
        INTwinsigma2        %[m^2], sigma^2 window interferers      (microphone spotformer)
        TARwinsigma2        %[m^2], sigma^2 window target           (microphone spotformer)      

        rebratio            %[],    ratio of direct to reverberant sound at mics
        NUMsigma2           %[m^2 (?)], variance of Rnum covariance matrix 
        Nsigma2             %[], variance of microphone self noise

        Ncontrol            %[-], the number of control points/person in the loudspeaker spotforming     
        
        flag_reverb         %[-], True if in simulated reverberant environment. False otherwise
        flag_real           %[-], True if in real environment. False otherwise
        flag_full_axis      %[-], True for full frequency axis [-Fs/2, Fs/2). False for [0, Fs/2]
        flag_nyquist        %[-], True if the nyquist bin for the microphone sptoformer should be set to zero artificially
        flag_BPF            %[-], True to apply a bandpass filter in the summation of the cost function of the loudspeaker spotformer
        flag_par_threshold  %[-], True to ignore low (<100 Hz) frequencies in the Par-function. I.e. the weights are set so high that they are effectively ignored
        flag_donut          %[-], True to integrate over a donut for the microphone array
    end

    methods
        function obj = Settings(c,  Nx, Ny, Nz, Nr, Ntheta, fs, window_length, pad_length, VAwinrad_x, VAwinrad_y, VAwinrad_z,...
                                VAwinrad_r, mu_r, INTwinrad, TARwinrad, rebratio, NUMsigma2, Nsigma2, Ncontrol, ...
                                flag_reverb, flag_real, flag_full_axis, flag_nyquist, flag_BPF, flag_par_threshold,...
                                flag_donut)
            if nargin == 0
                c = 342;
                Nx = 10;
                Ny = 10;
                Nz = 10;
                Nr = 10;
                Ntheta = 10;
                fs = 8000;
                window_length = 0.016;
                pad_length = 0.016;
                VAwinrad_x = 0.4;
                VAwinrad_y = 0.4;
                VAwinrad_z = 0.2;
                VAwinrad_r = 0.09;
                mu_r = 0.09;
                INTwinrad = 0.2;
                TARwinrad = 0.5;
                rebratio = 0;
                NUMsigma2 = 1e-10;
                Nsigma2 = 0;
                Ncontrol = 6;
                flag_reverb = false;        
                flag_real = false;
                      
                flag_full_axis = false;      
                flag_nyquist = false;       
                flag_BPF = false;            
                flag_par_threshold = false;
                flag_donut = true;
            end
            obj.c = c;
            obj.Nx = Nx;
            obj.Ny = Ny;
            obj.Nz = Nz;
            obj.Nr = Nr;
            obj.Ntheta = Ntheta;
            obj.fs = fs;
            
            obj.N_t = 2^nextpow2(floor(fs.*window_length));    
            obj.N_pad = 2^nextpow2(floor(fs.*pad_length));   
            obj.window_length_act = obj.N_t/fs;
            obj.pad_length_act = obj.N_pad/fs;
            
            obj.VAwinsigma2x = (VAwinrad_x/3)^2; 
            obj.VAwinsigma2y = (VAwinrad_y/3)^2;       
            obj.VAwinsigma2z = (VAwinrad_z/3)^2;       
            obj.VAwinsigma2r = (VAwinrad_r/3)^2;
            obj.mu_r = mu_r;
            obj.INTwinsigma2 = (INTwinrad/3)^2;         
            obj.TARwinsigma2 = (TARwinrad/3)^2;         

            obj.rebratio = rebratio;
            obj.NUMsigma2 = NUMsigma2;
            obj.Nsigma2 = Nsigma2;

            obj.Ncontrol = Ncontrol;         

            obj.flag_reverb = flag_reverb;    
            obj.flag_real = flag_real;       
            obj.flag_full_axis = flag_full_axis;      
            obj.flag_nyquist = flag_nyquist;       
            obj.flag_BPF = flag_BPF;             
            obj.flag_par_threshold = flag_par_threshold;
            obj.flag_donut = flag_donut;
        end
        function y = readAudio(obj, audiopath)
            [y, Fsamp] = audioread(audiopath);
            if  Fsamp ~= obj.fs
                disp("target sampling frequency (" + num2str(obj.fs) + ") and actual sampling frequency (" + num2str(Fsamp) + ") are unequal. Resampling...")
                y = resample(y, obj.fs, Fsamp);
            end
        end
    end
end
