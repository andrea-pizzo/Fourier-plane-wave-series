Fourier Plane-Wave Series Expansion for Holographic MIMO Communications
==================

This is a code package is related to the following article:

A. Pizzo, T. L. Marzetta and L. Sanguinetti, "Fourier Plane-Wave Series Expansion for Holographic MIMO Communications," in IEEE Transactions on Wireless Communications, vol. 21, no. 9, pp. 6890-6905, Sept. 2022. doi: 10.1109/TWC.2022.3152965.

The package contains a MATLAB script and related functions to compute the variances of Fourier random coefficients as described in Eq. (30) of the foregoing article. The computation is based on the power angle spectrum outlined in Eq. (31) and reproduces Figures 3 and 4.

The computed variances are essential for generating the Kronecker channel specified in Eq. (50) within Corollary 2. In a broader context, this script facilitates the generation of non-separable channels in accordance with Eq. (46) within Lemma 1. These non-separable channels would require a 6D power angle spectrum that depends jointly on the transmit-receive directions.

## Abstract of Article

Imagine a MIMO communication system that fully exploits the propagation characteristics offered by an electromagnetic channel and ultimately approaches the limits imposed by wireless communications. This is the concept of Holographic MIMO communications. Accurate and tractable channel modeling is critical to understanding its full potential. Classical stochastic models used by communications theorists are derived under the electromagnetic far-field assumption, i.e. planar wave approximation over the array. However, such assumption breaks down when electromagnetically large (compared to the wavelength) antenna arrays are considered. 

In this paper, we start from the first principles of wave propagation and provide a Fourier plane-wave series expansion of the channel response, which fully captures the essence of electromagnetic propagation in arbitrary scattering and is also valid in the (radiative) near-field. The expansion is based on the Fourier spectral representation and has an intuitive physical interpretation, as it statistically describes the angular coupling between source and receiver. When discretized uniformly, it leads to a low-rank semi-unitarily equivalent approximation of the electromagnetic channel in the angular domain. The developed channel model is used to compute the ergodic capacity of a point-to-point Holographic MIMO system with different degrees of channel state information.

## Content of Code Package

The package contains one Matlab script "main.m" and four additional Matlab functions "function_channelPAS.m", "function_channelVAR.m", "function_channelVAR_element.m", and "function_computeInt_nonIso.m", which are called by the main script.

See each file for further documentation. 


## License and Referencing

This code package is distributed under the GPLv2.0 license. If you utilize this code for research leading to publications, kindly acknowledge our original article mentioned above in your citations.
