function [var_channel_lm] = function_channelVAR_element(Lx,l,Ly,m,pas_channel)

%Compute the channel variances for every spherical surface element of the upper hemisphere.  
%%% INPUT VARIABLES
%Lx: array aperture along the x-axis in number of wavelengths. It must be an integer.
%l: index of the spherical surface element along the x-dimension.
%Ly: array aperture along the y-axis in number of wavelengths. It must be an integer.
%m: index of the spherical surface element along the y-dimension.
%pas_channel: function handle modeling the power angle spectrum of the channel. 
% It must subtend a unit area when integrated over the upper hemisphere.
%%% OUTPUT VARIABLES
%var_channel_lm: eigenvalue of the continuous channel computed over the (l,m)th 
% spherical surface element.

%Compute auxiliary parameters
a = l/Lx;
b = (l+1)/Lx;
c = m/Ly;
d = (m+1)/Ly;

%Compute the integral for any other regions
% First orthant %
if l>=0 && m>=0
    %Parametrize the azimuth angle
    phi1 = atan2(c,b);
    phi2 = min(atan2(c,a),atan2(d,b));
    phi3 = max(atan2(c,a),atan2(d,b));
    phi4 = atan2(d,a);
    %Compute the elevation angle integration interval
    theta_min_1 = @(p) asin(min(1,c./sin(p)));
    theta_max_1 = @(p) asin(min(1,b./cos(p)));
    %
    if l>=m
        theta_min_2 = @(p) asin(min(1,a./cos(p)));
        theta_max_2 = @(p) asin(min(1,b./cos(p)));
    else
        theta_min_2 = @(p) asin(min(1,c./sin(p)));
        theta_max_2 = @(p) asin(min(1,d./sin(p)));
    end
    %
    theta_min_3 = @(p) asin(min(1,a./cos(p)));
    theta_max_3 = @(p) asin(min(1,d./sin(p)));
    
    % Second orthant %
elseif l<0 && m>=0
    %Parametrize the azimuth angle
    phi1 = pi - atan2(d,abs(b));
    phi2 = min(pi-atan2(c,abs(b)),pi-atan2(d,abs(a)));
    phi3 = max(pi-atan2(c,abs(b)),pi-atan2(d,abs(a)));
    phi4 = pi - atan2(c,abs(a));
    %Compute the elevation angle integration interval
    theta_min_1 = @(p) asin(min(1,abs(b)./cos(pi-p)));
    theta_max_1 = @(p) asin(min(1,d./sin(pi-p)));
    %
    if abs(l)>m
        theta_min_2 = @(p) asin(min(1,abs(b)./cos(pi-p)));
        theta_max_2 = @(p) asin(min(1,abs(a)./cos(pi-p)));
    elseif abs(l)<=m
        theta_min_2 = @(p) asin(min(1,abs(c)./sin(pi-p)));
        theta_max_2 = @(p) asin(min(1,abs(d)./sin(pi-p)));
    end
    %
    theta_min_3 = @(p) asin(min(1,c./sin(pi-p)));
    theta_max_3 = @(p) asin(min(1,abs(a)./cos(pi-p)));
    
    % Third orthant %
elseif l<0 && m<0
    %Parametrize the azimuth angle
    phi1 = pi + atan2(abs(d),abs(a));
    phi2 = min(pi+atan2(abs(d),abs(b)),pi+atan2(abs(c),abs(a)));
    phi3 = max(pi+atan2(abs(d),abs(b)),pi+atan2(abs(c),abs(a)));
    phi4 = pi + atan2(abs(c),abs(b));
    %Compute the elevation angle integration interval
    theta_min_1 = @(p) asin(min(1,abs(d)./sin(p-pi)));
    theta_max_1 = @(p) asin(min(1,abs(a)./cos(p-pi)));
    %
    if abs(l)>=abs(m)
        theta_min_2 = @(p) asin(min(1,abs(b)./cos(p-pi)));
        theta_max_2 = @(p) asin(min(1,abs(a)./cos(p-pi)));
    else
        theta_min_2 = @(p) asin(min(1,abs(d)./sin(p-pi)));
        theta_max_2 = @(p) asin(min(1,abs(c)./sin(p-pi)));
    end
    %
    theta_min_3 = @(p) asin(min(1,abs(b)./cos(p-pi)));
    theta_max_3 = @(p) asin(min(1,abs(c)./sin(p-pi)));
    
    % Fourth orthant %
elseif l>=0 && m<0
    %Parametrize the azimuth angle
    phi1 = 2*pi - atan2(abs(c),a);
    phi2 = min(2*pi-atan2(abs(c),b),2*pi-atan2(abs(d),a));
    phi3 = max(2*pi-atan2(abs(c),b),2*pi-atan2(abs(d),a));
    phi4 = 2*pi - atan2(abs(d),b);
    if d==0 && a==0 %undetermined form inside atan2
        phi2 = 3/2*pi;
        phi3 = 2*pi - atan2(abs(c),b);
        phi4 = 2*pi;
    end
    %Compute the elevation angle integration interval
    theta_min_1 = @(p) asin(min(1,a./cos(2*pi-p)));
    theta_max_1 = @(p) asin(min(1,abs(c)./sin(2*pi-p)));
    %
    if l>=abs(m)
        theta_min_2 = @(p) asin(min(1,a./cos(2*pi-p)));
        theta_max_2 = @(p) asin(min(1,b./cos(2*pi-p)));
    else
        theta_min_2 = @(p) asin(min(1,abs(d)./sin(2*pi-p)));
        theta_max_2 = @(p) asin(min(1,abs(c)./sin(2*pi-p)));
    end
    %
    theta_min_3 = @(p) asin(min(1,abs(d)./sin(2*pi-p)));
    theta_max_3 = @(p) asin(min(1,b./cos(2*pi-p)));
    %
    if l==0 && m==-1
        theta_min_2 = @(p) zeros(size(p));
        %
        theta_min_3 = @(p) zeros(size(p));
    end
    
end

%Compute the double integrals
var_channel_lm_1 = function_computeInt_nonIso(theta_min_1,theta_max_1,phi1,phi2,pas_channel);
var_channel_lm_2 = function_computeInt_nonIso(theta_min_2,theta_max_2,phi2,phi3,pas_channel);
var_channel_lm_3 = function_computeInt_nonIso(theta_min_3,theta_max_3,phi3,phi4,pas_channel);

%Sum up the three contributions
var_channel_lm = var_channel_lm_1 + var_channel_lm_2 + var_channel_lm_3;

end

