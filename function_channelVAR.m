function [var_channel] = function_channelVAR(Lx,Ly,pas_channel)

%Generate the channel variances for a given power angle spectrum.
%
%%% INPUT VARIABLES
%Lx: array aperture in number of wavelengths along the x-axis. It must be an integer.
%Ly: array aperture in number of wavelengths along the y-axis. It must be an integer.
%pas_channel: power angle spectrum function of the channel. Its integration over the unit upper 
% hemisphere must yield one.
%
%%% OUTPUT VARIABLES
%var_channel: (2*Lx x 2*Ly) matrix with nonnegative entries modeling the approximate 
% eigenvalues of the continuous channel. The sum of the matrix entries must yield 1.

%Spatial DOF per dimension (rectangular embedding)
ny = 2*Ly;
nx = 2*Lx;

%Compute the channel variances
if Ly==0 %ULA
    lx_vec = (-Lx:1:Lx-1); %spatial frequencies
    ly_vec = (-1:1:0)';
    var_channel = NaN*ones(2,nx);
    for ll=1:nx
        l = lx_vec(ll);
        for mm=1:2
            m = ly_vec(mm);
            var_channel(end-mm+1,ll) = function_channelVAR_element(Lx,l,1,m,pas_channel);
        end
    end
    var_channel = sum(var_channel,1); %compute variances over the x-dimension only
else %UPA
    lx_vec = (-Lx:1:Lx-1); 
    ly_vec = (-Ly:1:Ly-1)';
    var_channel = NaN*ones(ny,nx);
    for ll=1:nx
        l = lx_vec(ll);
        for mm=1:ny
            m = ly_vec(mm);
            var_channel(end-mm+1,ll) = function_channelVAR_element(Lx,l,Ly,m,pas_channel);
        end
    end
end

%normalize the variances to unit power
var_channel(var_channel<0) = 0;
var_channel = var_channel/sum(var_channel(:));

end