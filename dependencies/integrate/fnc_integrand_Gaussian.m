%Author:    Dimme de Groot
%Date:      May 2024
%Descr:     This function is used for computing the spatial covariance matrices. 
%               out = fnc_integrand_Gaussian(x_s1, x_s2, x_r, k, mu_r, sigma_x, sigma_y, sigma_z, flag_spatial_weight)
%           If flag_spatial_weight==false     
%               out = 1/(sqrt{(2*pi)^3}*sigma_x*sigma_y*sigma_z) * g(x_s1, x_r, k)*conj(g(x_s2, x_r, k))
%           If flag_spatial_weight==true, a Gaussian is added:
%               out = 1/(sqrt{(2*pi)^3}*sigma_x*sigma_y*sigma_z) * g(x_s1, x_r, k)*conj(g(x_s2, x_r, k))*exp(-0.5 {[(x_r(1)-mu(1))/sigma_x]^2+[(x_r(2)-mu(2))/sigma_y]^2+[(x_r(3)-mu(3))/sigma_z]^2})
%
%           Here, g is the greens function solution to the wave equation
%                      
%           Inputs: 
%               - x_s1, x_s2 the source locations. 3x1 real vectors
%               - x_r the receiver location. 3xN real vector --> the output is computed for each N
%               - mu_r the mean of the Guassian. 3x1 real vector
%               - sigma_x, sigma_y, sigma_z the standard deviations. Real scalar
%               - k the wave number. Real scalar
%               - flag_spatial_weight selects if the Gaussian weighting should be added. Note that, if this is false, the normal;isation term is still(!!) added!
%           Output:
%               - Out: 1xN complex vector. Each element out(1,i) is the result evaluated for receiver location x_r(:,i);

function out = fnc_integrand_Gaussian(x_s1, x_s2, x_r, k, mu_r, sigma_x, sigma_y, sigma_z, flag_spatial_weight)
    %The inner term of the integrand
    out = fnc_green_wavefunction(x_s1, x_r, k).*conj(fnc_green_wavefunction(x_s2, x_r, k));
    
    %Normalisation term of (Gaussian) weighting
    normConst = sqrt((2*pi)^3)*sigma_x*sigma_y*sigma_z;
    
    %Gaussian (spatial) weighting
    if flag_spatial_weight     %Add the spatial weighting term
        diff = x_r - mu_r;
        out = out.*exp(-0.5*( (diff(1,:)/sigma_x).^2 + (diff(2,:)/sigma_y).^2  + (diff(3,:)/sigma_z).^2 ));
    end

    %Output should be normalised by normalisation term
    out = out/normConst;
end

