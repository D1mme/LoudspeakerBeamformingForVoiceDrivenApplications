%Author:    Dimme de Groot
%Date:      March 2024
%Descr:     This class is used to describe a room. 
%           Each room has:
%               a speed of sound c              [m/s]  
%               a speaker location S            [m],    (S of size Ns x 3)
%               a microphone location R         [m],    (R of size Nr x 3)
%               a person location P             [m],    (P of size 1 x 3)
%               mean microphone location Rbar   [m],    (Rbar of size 1 x 3)
%               walls of length L               [m],    (L = [Lx, Ly, Lz])   
%               reflection coefficients beta    [-],    ([bx, bx, by, by, bz, bz])
%               number of speakers Ns           [-],        
%               number of microphones Nr        [-],   
%               Number of persons Np            [-],
%               An estimated T60 time           [s],    see [1].
%Functions: 
%   room = Room(R, S, L, beta, sound_velocity): this creates a class instant room, where sound velocity can be excluded
%   room.plotRoom(top_view_flag): this function provides a basic plot of the room. Set top_view_flag = true for a top view (x,y plane).             
%   [hLoudspeaker, hPerson] = room.comp_transfer(fs, mean_flag), with fs the sampling frequency and mean_flag is false to get the trasnfer to R and true to get the transfer to Rbar    
%       hLoudspeaker of size (fs*T60, Ns, {1 or Nr})
%       hPerson of size (fs*T60, Np, {1 or Np})
%
%Dependencies:  
%   using room.comp_transfer requires rir_generator.mex, see [1].           
%
%Sources:
%   [1] E. Habets, RIR generator, https://github.com/ehabets/RIR-Generator

classdef Room
    properties
        L           %[m], wall lengths (x, y, z)
        beta        %[-], reflection coefficients    
        R           %[m], microphone locations
        S           %[m], loudspeaker locations
        P           %[m], person location (person giving voice commands)

        Rbar        %[m], mean mic location
        NN          %[-], the microphone index corresponding to the nearest neighbouring microphone of the person
        Ns          %[-], number of loudspeakers
        Nr          %[-], number of mics
        Np          %[-], number of persons

        c           %[m/s], speed of sound
        T60         %[s], the "T60" time of the room. Computed using estimate described in [1]      
    end

    methods
        function obj = Room(R, S, P, L, beta, sound_velocity)
            if nargin==4
                obj.c = 342;
            else
                obj.c = sound_velocity;  
            end
            obj.R = R;
            obj.S = S;
            obj.P = P;  
            obj.L = L;
            obj.beta = beta;
            obj.Ns = size(obj.S,1);
            obj.Nr = size(obj.R,1);
            obj.Np = size(obj.P,1);
            obj.T60 = methodGetT60(obj);
            obj.Rbar = methodGetRbar(obj);
            [~, obj.NN] = min(vecnorm(obj.R-obj.P,2,2));
        end 
        function T60 = methodGetT60(obj)    %see [1]
            aS_x = obj.L(2)*obj.L(3)*(1-obj.beta(1)^2+1-obj.beta(2)^2);
            aS_y = obj.L(1)*obj.L(3)*(1-obj.beta(3)^2+1-obj.beta(4)^2);
            aS_z = obj.L(1)*obj.L(2)*(1-obj.beta(5)^2+1-obj.beta(6)^2);
            V = obj.L(1)*obj.L(2)*obj.L(3);
            T60 = 24*log(10)*V/(obj.c*(aS_x+aS_y+aS_z));
        end
        function Rbar = methodGetRbar(obj)
            Rbar = mean(obj.R); %mean gets taken per column
        end
        function plotRoom(obj, top_view_flag, annotate_flag)
            %plotcolors: https://personal.sron.nl/~pault/
            %high_contrast  = [[221 170 51]; [187 85 102]; [0 68 136]]/255; 
            if nargin==1
                top_view_flag = false;
                annotate_flag = false;
            end
            if nargin==2
                annotate_flag = false;
            end
            if top_view_flag
                figure
                %draw sources
                scatter(obj.R(:,1), obj.R(:,2), "filled")
                hold on
                scatter(obj.S(:,1), obj.S(:,2), "filled")
                scatter(obj.P(:,1), obj.P(:,2), "filled")
                if annotate_flag
                    text(obj.R(:,1), obj.R(:,2), num2str((1:obj.Nr)'));
                    text(obj.S(:,1), obj.S(:,2), num2str((1:obj.Ns)'));
                end
                %draw walls room
                plot([0 0], [0, obj.L(2)], 'black',  'linewidth',2)
                plot([0, obj.L(1)], [0, 0], 'black', 'linewidth',2)
                plot([obj.L(1) obj.L(1)], [0, obj.L(2)], 'black',  'linewidth',2)
                plot([0, obj.L(1)], [obj.L(2), obj.L(2)], 'black', 'linewidth',2)
                
                xlim([0-0.3 obj.L(1)+0.3])
                ylim([0-0.3, obj.L(2)+0.3])
            else
                high_contrast  = [[221 170 51]; [187 85 102]; [0 68 136]]/255; 

                figure
                %draw mics and loudspeakers
                scatter3(obj.R(:,1), obj.R(:,2), obj.R(:,3), "filled", 'MarkerFaceColor', high_contrast(1,:), 'MarkerEdgeColor', high_contrast(1,:))
                hold on
                scatter3(obj.S(:,1), obj.S(:,2), obj.S(:,3), "filled", 'MarkerFaceColor', high_contrast(2,:), 'MarkerEdgeColor', high_contrast(2,:)) 
                scatter3(obj.P(:,1), obj.P(:,2), obj.P(:,3), "filled", 'MarkerFaceColor', high_contrast(3,:), 'MarkerEdgeColor', high_contrast(3,:))
             
                if annotate_flag
                    text(obj.R(:,1), obj.R(:,2), obj.R(:,3)+0.1, num2str((1:obj.Nr)'));
                    text(obj.S(:,1), obj.S(:,2), obj.S(:,3)+0.1, num2str((1:obj.Ns)'));
                end
                
                %draw some walls room
                plot3([0, 0], [obj.L(2), obj.L(2)], [0, obj.L(3)],  'black',  'linewidth',2)
                plot3([0 0], [0, obj.L(2)], [0,0], 'black',  'linewidth',2)
                plot3([0, obj.L(1)], [0, 0], [0,0], 'black', 'linewidth',2)

                xlim([0 obj.L(1)])
                ylim([0, obj.L(2)])
                zlim([0, obj.L(3)])
            end
            grid on
            axis equal
            legend('Microphone', 'Loudspeaker', 'Person','Location','northeast')
            xlabel('xr [m]')
            ylabel('yr [m]')
            zlabel('zr [m]')
            set(gcf,'renderer','Painters')
        end
        function [hLoud, hPerson, hLoudPerson] = comp_transfer(obj, fs, mean_flag)
            if mean_flag
                Rmic = obj.Rbar;
            else
                Rmic = obj.R;
            end
            mtype = 'omnidirectional';  % Type of microphone
            order = -1;                 % Reflection order
            dim = 3;                    % Room dimension
            orientation = 0;            % Microphone orientation (rad)
            hp_filter = 0;              % Enable high-pass filter
        
            Nh = ceil(obj.T60*fs);
            Nmic = size(Rmic,1); 
            hLoud = zeros(Nh, obj.Ns, Nmic);
            hPerson = zeros(Nh, obj.Np, Nmic);
            hLoudPerson = zeros(Nh, obj.Ns, obj.Np);
            for i = 1:obj.Ns
                for j=1:Nmic
                    evalc('hLoud(:,i,j) =rir_generator(obj.c, fs, Rmic(j,:), obj.S(i,:), obj.L, obj.beta, Nh, mtype, order, dim, orientation, hp_filter);');
                    %hLoud(:,i,j) = rir_generator_x(obj.c, fs, Rmic(j,:), obj.S(i,:), obj.L, obj.beta, Nh);
                end
            end

            for i=1:obj.Np
                for j=1:Nmic
                    evalc('hPerson(:,i,j) = rir_generator(obj.c, fs, Rmic(j,:), obj.P(i,:), obj.L, obj.beta, Nh, mtype, order, dim, orientation, hp_filter);');
                    %hPerson(:,i,j) = rir_generator_x(obj.c, fs, Rmic(j,:), obj.P(i,:), obj.L, obj.beta, Nh);
                end
            end

            for i=1:obj.Np
                for j=1:obj.Ns
                    evalc('hLoudPerson(:,j,i) = rir_generator(obj.c, fs, obj.S(j,:), obj.P(i,:), obj.L, obj.beta, Nh, mtype, order, dim, orientation, hp_filter);');
                    %hPerson(:,i,j) = rir_generator_x(obj.c, fs, Rmic(j,:), obj.P(i,:), obj.L, obj.beta, Nh);
                end
            end
        end
    end
end
