# TWC23
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
