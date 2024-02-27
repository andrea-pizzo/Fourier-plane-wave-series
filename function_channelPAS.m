function [pas_channel] = function_channelPAS(circular_var,mean_rad_theta,mean_rad_phi)

%Generate the power angle spectrum of a channel according to the mixture of Von-Mises Fisher
% distributions. 
%
%%% INPUT VARIABLES (N=number of clusters)
%circular_var: (N x 1) specifies the angular concentration around the mean direction. It must be 
% within (0,1] with 1 corresponding to isotropic scattering and 0 to an impulsive scattering.
%mean_rad_theta: (N x 1) specifies the elevation angle of the mean direction wrt the z-axis. 
% It must be within [0,pi/2) for upgoing propagation.
%mean_rad_phi: (N x 1) specifies the azimuth angle of the mean direction wrt the x-axis. 
% It must be within [0,2pi).
%
%%% OUTPUT VARIABLES
%pas_channel: power angle spectrum function of the channel. Its integration over the unit upper 
% hemisphere must yield one.

%Number of clusters
numOfClusters = length(circular_var);

%initialize
pas_channel = @(p,t) 0;
for n=1:numOfClusters

    circular_var_n = circular_var(n);
    mean_rad_theta_n = mean_rad_theta(n);
    mean_rad_phi_n = mean_rad_phi(n);

    %Generate the angular density function for the nth cluster
    if circular_var_n==1 %isotropic scattering
        alphaval_n = 0;
        pas_channel_n = @(p,t) 1/(2*pi);
    elseif circular_var_n==0 %impulsive scattering
        alphaval_n = infty;
        pdf_n = NaN;
    else
        %Compute the auxiliary concentration parameter
        myfun = @(a) 1-(coth(a)-1/a)^2-circular_var_n; %v=1-rho Eq.(3.4.11), rho=coth(a)-1/a Eq.(9.3.9) - Prof. Mardia book
        a0 = [1e-6, 1e6]; %search interval restricted to positive numbers only
        [alphaval_n fval exitflag output] = fzero(myfun,a0);
        %Generate the angular density
        cos_gamma_fun = @(p,t) sin(t)*sin(mean_rad_theta_n).*cos(p-mean_rad_phi_n) + cos(t)*cos(mean_rad_theta_n);
        pas_channel_n = @(p,t) alphaval_n*exp(alphaval_n*cos_gamma_fun(p,t))/(4*pi*sinh(alphaval_n));
    end

    %Generate the mixture of density functions
    pas_channel = @(p,t) pas_channel(p,t) + pas_channel_n(p,t)/numOfClusters;

end

%Multiply the angular density by the spherical Jacobian
pas_channel = @(p,t) pas_channel(p,t).*sin(t);

%Check that the density function subtends unit area in the upper hemisphere
pas_channel_normalization = integral2(pas_channel,0,2*pi,0,pi/2)

end