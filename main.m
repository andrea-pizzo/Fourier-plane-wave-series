% MATLAB Script
%
% This MATLAB script is designed to compute the variances of Fourier
% random coefficients as described in Eq. (30). The computation is based on
% the power angle spectrum outlined in Eq. (31), as presented in the article:
%
%A. Pizzo, T. L. Marzetta and L. Sanguinetti, "Fourier Plane-Wave Series Expansion for 
% Holographic MIMO Communications," in IEEE Transactions on Wireless Communications, 
% vol. 21, no. 9, pp. 6890-6905, Sept. 2022. doi: 10.1109/TWC.2022.3152965.
%
% The computed variances are essential for generating the Kronecker channel
% specified in Eq. (50) within Corollary 2. In a broader context, this
% script facilitates the generation of non-separable channels in accordance
% with Eq. (46) within Lemma 1. These non-separable channels require a
% 6D power angle spectrum that depends jointly on the transmit-receive directions. 
%
% License: This code is distributed under the GPLv2.0 license. If you utilize
% this code for research leading to publications, kindly acknowledge our
% original article mentioned above in your citations.

%Empty workspace and close figures
clear
close all
clc

%%%% Parameters %%%% 
%Array aperture (in number of wavelengths)
Lx = 10; %integer number
Ly = 10; %integer number
%Power angle spectrum (3D Von Mises-Fisher)
circular_var_vec = [.01,.005]; %circular variance \in(0,1] (1=isotropic, regardless of mean direction)
mean_theta_deg_vec = [30,10]; % mean direction elevation (in degrees) \in[0,90)
mean_phi_deg_vec = [15,180];  % mean direction azimuth (in degrees) \in[0,360)

%%%% Simulation %%%%
%2D spatial frequencies
lx_vec = (-Lx:1:Lx-1);
ly_vec = (-Ly:1:Ly-1);

%Mean direction (in radians)
mean_theta_vec = mean_theta_deg_vec/180*pi;
mean_phi_vec = mean_phi_deg_vec/180*pi;

%Generate the power angle spectrum
pas_channel = function_channelPAS(circular_var_vec,mean_theta_vec,mean_phi_vec);

%Compute the variances of Fourier coefficients
var_channel = function_channelVAR(Lx,Ly,pas_channel);

%plot the power angle spectrum in Figure 4
theta_vec = linspace(0,pi/2,100);
phi_vec = linspace(0,2*pi,100);
[theta,phi] = meshgrid(theta_vec,phi_vec);
[x,y,z] = sph2cart(phi,pi/2-theta,1); %NOTE: elevation angle is defined wrt the z-axis 
figure;
surf(x,y,z,10*log10(pas_channel(phi,theta))); hold on;
caxis([-30 0]);
colorbar; shading interp; daspect([1 1 1]); axis tight;
xlabel('$y$','Interpreter','Latex');
ylabel('$x$','Interpreter','Latex');
zlabel('$z$','Interpreter','Latex');
title('$A^2_{\rm R}(\theta,\phi)$','Interpreter','Latex')
grid on; box on;
view([115 40]);
set(gca,'FontSize',20);

%plot the channel variances (dB) in Figure 4
figure;
[X,Y] = meshgrid(lx_vec,ly_vec);
surf(X,Y,10*log10(var_channel/max(var_channel(:))));
caxis([-60 0]);
colorbar
xlabel('$\ell_x$','Interpreter','Latex');
ylabel('$\ell_y$','Interpreter','Latex');
xlim([-Lx Lx-1])
ylim([-Ly Ly-1])
zlabel('$\sigma^2_{\rm R}(\ell_x,\ell_y)$ (dB)','Interpreter','Latex');
grid on; box on;
view([0 90]);
set(gca,'FontSize',20);

