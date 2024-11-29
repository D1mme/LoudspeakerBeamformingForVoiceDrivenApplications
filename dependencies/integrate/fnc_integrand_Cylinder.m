%Author:    Dimme de Groot
%Date:      May 2024
%Descr:     This function is used for computing the spatial covariance matrices. 
%               out = fnc_integrand_Cylinder(x_s1, x_s2, x_bar, r, theta, z, k, mu_r, mu_z, sigma_r, sigma_z, flag_spatial_weight)                     
%           Inputs: 
%               - x_s1, x_s2 the source locations. 3x1 real vectors
%               - x_bar the receiver location. 3x1 real vector. Effectively, this represent the coordinate around which the cylinder is centered 
%		        - r, theta, z: 1xN real vectors. These vectors specify the coordinate of the receiver when centered around x_bar (i.e z=0 effectively represents z'=x_bar(3))
%               - k the wave number. Real scalar
%               - mu_r the mean value of the radius
%               - sigma_r, sigma_z: the standard deviation in the weighting of the coordinates r and z
%               - flag_spatial_weight selects if the Gaussian weighting terms be added. Note that, if this is false, the normal;isation term is still(!!) added!
%           Output:
%               - Out: 1xN complex vector. Each element out(1,i) is the result evaluated for receiver location [r*cos(i), r*sin(i), z(i)]+Rbar

function out = fnc_integrand_Cylinder(x_s1, x_s2, x_bar, r, theta, z, k, mu_r, mu_z, sigma_r, sigma_z, flag_spatial_weight)
    if mu_z ~= 0
        disp("fnc_integrand_Cylinder: you might want to check mu_z");
    end

    %The inner term of the integrand
    x_r(1,:) = r.*cos(theta)+x_bar(1);
    x_r(2,:) = r.*sin(theta)+x_bar(2);
    x_r(3,:) = z + x_bar(3);
    out = fnc_green_wavefunction(x_s1, x_r, k).*conj(fnc_green_wavefunction(x_s2, x_r, k)).*r;

    %Normalisation term of (Gaussian) weighting
    normConst = 2*pi*sigma_r*sigma_z*2*pi; %Two from the Gaussians, one from the uniform

    %Spatial weighting
    if flag_spatial_weight     %Add the spatial weighting term
        diff = r - mu_r;
        out = out.*exp(-0.5*(z/sigma_z).^2).*exp(-0.5*(diff/sigma_r).^2).*fnc_step(r);
    end

    %Output should be normalised by normalisation term
    out = out/normConst;

    function u = fnc_step(r)
        u = ones(size(r));
        u(r<0) = 0;
    end
end

