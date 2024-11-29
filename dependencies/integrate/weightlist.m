%Author:    Dimme de Groot
%Date:      May 2024
%Descr:     Inefficient but easy code to generate a list of coordinates and the corresponding weights. Used for numerical integration
%           Inputs:
%               w_x, coor_x, w_y, coor_y, w_z, coor_z
%                   w_{x,y,z} are weightvectors of length N_{x,y,z}. Orientation does not matter
%                   coor_{x,y,z} are the corresponding coordinates of length  N_{x,y,z}. Orientation does not matter
%           Outputs:
%               weight_list, coor_list
%                   coor_list is a list of coordinates of size [D, N_x*N_y*N_z]. Here D=1 if you only gave x, D=2 for x and y, D=3 for x and y and z
%                   weight_list is a list of weights of size [1, N_x*N_y*N_z] with the weight of each coordinate of the coordinate list. 

function [weight_list, coor_list] = weightlist(w_x, coor_x, w_y, coor_y, w_z, coor_z, verbose)
    if nargin < 7
        verbose = false;
    end
    if nargin == 0 || nargin == 1 || nargin == 3 || nargin == 5 || nargin > 7
        disp('weightlist: invalid number of inputs');
    elseif nargin == 2
        if verbose
            disp('weightlist: detected 1 dimensional coordinate system')
        end
        D = 1;
    elseif nargin == 4
        if verbose
            disp('weightlist: detected 2 dimensional coordinate system')
        end       
        D = 2;
        coor_z = [];
        w_z = [];
    elseif nargin == 6
        if verbose
            disp('weightlist: detected 3 dimensional coordinate system')
        end
        D = 3;
    end
    
    %Define lengths and placeholds
    if D == 1
        Nx = length(w_x);
        Ny = 1;
        Nz = 1;
        coor_list = zeros(D,Nx);
        weight_list = zeros(1,Nx);
    elseif D==2
        Nx = length(w_x);
        Ny = length(w_y);
        Nz = 1;
        coor_list = zeros(D,Nx*Ny);
        weight_list = zeros(1,Nx*Ny);
    elseif D==3
        Nx = length(w_x);
        Ny = length(w_y);
        Nz = length(w_z);
        coor_list = zeros(D,Nx*Ny*Nz);
        weight_list = zeros(1,Nx*Ny*Nz);
    end
    
    %Very inefficent but bugsafe: loop through coordinates
    count = 1;
    for i=1:Nx
        for j=1:Ny
            for k=1:Nz
                if D==1
                    coor_list(:,count) = coor_x(i);
                    weight_list(1,count) = w_x(i);
                elseif D==2
                    coor_list(:,count) = [coor_x(i); coor_y(j)];
                    weight_list(1,count) = w_x(i)*w_y(j);
                elseif D==3
                    coor_list(:,count) = [coor_x(i); coor_y(j); coor_z(k)];
                    weight_list(1,count) = w_x(i)*w_y(j)*w_z(k);
                end
                count = count+1;
            end
        end
    end
end