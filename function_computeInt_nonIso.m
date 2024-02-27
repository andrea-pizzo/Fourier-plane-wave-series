function [int] = function_computeInt_nonIso(theta_min,theta_max,phi_min,phi_max,pas_channel)

%Compute the channel variances for every spherical surface element of the upper hemisphere.
%
%%% INPUT VARIABLES
%theta_min: lower limit integration region in terms of the elevation angle.
%theta_max: upper limit integration region in terms of the elevation angle.
%phi_min: lower limit integration region in terms of the azimuth angle.
%phi_max: upper limit integration region in terms of the azimuth angle.
%pas_channel: function handle modeling the power angle spectrum of the channel.
% It must subtend a unit area when integrated over the upper hemisphere.
%
%%% OUTPUT VARIABLES
%int: auxiliary integral needed to compute the channel variances wrt the (l,m)th
% spherical surface element.

int = real(integral2(pas_channel,phi_min,phi_max,theta_min,theta_max,'method','iterated'));

end