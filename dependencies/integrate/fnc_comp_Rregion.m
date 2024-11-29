%Author:    Dimme de Groot
%Date:      March 2024
%Descr:     Computes the covariance matrices over the regions. (Excluding the isotropic and numerical covariance matrices)     

function R_region = fnc_comp_Rregion(N, k_ax, x, mu, sigma2, cylinder_flag, flag_clenshaw_curtis)
    if nargin == 5
        cylinder_flag = false;
        flag_clenshaw_curtis = true;
    end
    if nargin == 6
        flag_clenshaw_curtis = true;
    end
    if isscalar(sigma2)
        %If this happens we are doing microphone spotforming
        sigma2 = [sigma2, sigma2, sigma2, sigma2];
        mu = [mu,0];
        cylinder_flag = false; 
        N = [N, N, N, N, N];
    end
    Ns = size(x,1);  
    Nk = length(k_ax);
    
    Nx = N(1); Ny = N(2); Nz = N(3); Nr = N(4); Ntheta = N(5);
    sigma = sqrt(sigma2);
    sigma_x = sigma(1); sigma_y = sigma(2); sigma_z = sigma(3); sigma_r = sigma(4);
    x_bar = mu(1:3); mu_r = mu(4);
    
    R_region = zeros(Ns, Ns, Nk);
    
    if cylinder_flag
        %Donutshaped weighting for the loudspeaker spotformer

        if ~flag_clenshaw_curtis
            disp("Hermquad is not supported due to licensing")
        else
            %Compute quadrature points
            [r_samp, w_r] = clenquad(Nr, max(0, mu_r-3.3*sigma_r), 3.3*sigma_r+mu_r); %From -3.3sigma --> 3.3sigma: the 3sigma radius and a bit. I thought -3.3sigma corresponds to -60 dB, but (iirc) it's about 3.7sigma.  
            [theta_samp, w_theta] = clenquad(Ntheta, 0, 2*pi);
            [z_samp, w_z] = clenquad(Nz, -3.3*sigma_z, 3.3*sigma_z);
            [weight_list, coor_list] = weightlist(w_r, r_samp, w_theta, theta_samp, w_z, z_samp);
        end

        %Compute estimate PSD matrix by making a weighted sum over the quadrature points. Performance of this part can likel 
        for i=1:Ns
            for j=1:Ns
                for k=1:Nk
                    kk = k_ax(k);          
                    if ~flag_clenshaw_curtis %Gauss-Hermite: not supported due to licensing
                        out = fnc_integrand_Cylinder(x(i,:)', x(j,:)', x_bar', coor_list(1,:), coor_list(2,:), coor_list(3,:), kk, mu_r, 0, sigma_r, sigma_z, false);
                        R_region(i,j,k) = sum(2*sigma_r*sigma_z*weight_list.*out);
                    else %Clenshaw-curtis 
                        out = fnc_integrand_Cylinder(x(i,:)', x(j,:)', x_bar', coor_list(1,:), coor_list(2,:), coor_list(3,:),kk, mu_r, 0, sigma_r, sigma_z, true);
                        R_region(i,j,k) = sum(weight_list.*out);
                    end 
                end
            end
        end
    else
        %Gaussian weighting for the loudspeaker spotformer. We always have Gaussian weighting for the microphone spotformer
        if ~flag_clenshaw_curtis
            disp("Hermquad is not supported due to licensing")
        else
            %Compute quadrature points
            [x_quad, w_x] = clenquad(Nx, -3.3*sigma_x+x_bar(1), 3.3*sigma_x+x_bar(1));
            [y_quad, w_y] = clenquad(Ny, -3.3*sigma_y+x_bar(2), 3.3*sigma_y+x_bar(2)); 
            [z_quad, w_z] = clenquad(Nz, -3.3*sigma_z+x_bar(3), 3.3*sigma_z+x_bar(3));
            [weight_list, coor_list] = weightlist(w_x, x_quad, w_y, y_quad, w_z, z_quad);
        end

        for i=1:Ns
            for j=1:Ns
                for k=1:Nk
                    kk = k_ax(k);          
                    if ~flag_clenshaw_curtis %Gauss-Hermite: not supported due to licensing
                        out = fnc_integrand_Gaussian(x(i,:)', x(j,:)', coor_list, kk, x_bar', sigma_x, sigma_y, sigma_z, false);
                        R_region(i,j,k) = sum(weight_list.*out)*(sqrt(2)*sigma_x)*(sqrt(2)*sigma_y)*(sqrt(2)*sigma_z);
                    else %Clenshaw-curtis 
                        out = fnc_integrand_Gaussian(x(i,:)', x(j,:)', coor_list, kk, x_bar', sigma_x, sigma_y, sigma_z, true);
                        R_region(i,j,k) = sum(weight_list.*out);   
                    end 
                end
            end
        end
    end
end
