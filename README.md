# Framework for Terahertz Wireless Systems
This repository provides a Matlab-based framework to perform numerical analysis on channel models obtained upon combining small-scale fading and antenna misalignment effects. This analysis can also be verified using Monte-Carlo Simulations. We have also included code to particularly work with the Fluctuating Two Ray (FTR distribution), which is a recent addition to the channel models of Terahartz Wireless Systems. To use this code, clone this repository to your project folder.

## Why Fluctuating Two Ray?
A recent measurement campaign was conducted to measure the THz channel with operation at 304.2 GHz and 8 GHz bandwidth. Gradient descent was used to fit the distribution parameters for Gaussian, Rician, Nakagami-m and FTR distributions. FTR distribution performed a much better fit than other models. Morever, FTR distribution is a generalization of a lot more commonly used probability distributions for small scale fading.

## Fundamentals of FTR:
The small-scale fluctuations in the amplitude of a signal transmitted over a wireless channel can be modeled by the superposition of a set of N dominant waves, referred to as specular components, to which other diffusely propagating waves are added.

![image](https://user-images.githubusercontent.com/65544914/117104407-b26cac00-ad99-11eb-94af-76c602f3af47.png)

X + jY is a complex Gaussian random variable, such that X and Y follow a Normal distribution with zero mean and standard deviation sigma. The random phase follows a Uniform distribution between zero and 2 pi. The amplitudes Vn are constants.

## FTR Probabability Distribution:

![image](https://user-images.githubusercontent.com/65544914/117105154-1348b400-ad9b-11eb-9bb4-012dc871accf.png)

## Details on the use of this code:

- [get_ftr_pdf.m](get_ftr_pdf.m) : Returns a function handle to the FTR PDF. Can be used to verify integrals and for numerical analysis.
- [get_combined_pdf.m](get_combined_pdf.m) :  Returns a function handle to the Combined FTR + Misalignment errors PDF.
- [d_calc.m](d_calc.m) : Helper function to calculate the value of d as given in the FTR Probability Distribution expression above.
- [my_asso_legendre_func.m](my_asso_legendre_func.m) : This function is be used to calculate Associated Legendre function for real parameters and complex argument. \[For details, refer to I. S. Gradshteyn and I. M. Ryzhik, Table of Integrals, Series, and Products, 7th ed. Academic Press, 2007.\]   
- [THz_Pathloss.m](THz_Pathloss.m) : Helper function to calculate the deterministic path loss due to molecular absorption.
- [randpdf.m](randpdf.m) : Helper function to generate pseudo-random numbers from a user defined distribution. \[Credits: Adam Nies³ony, Opole University of Technology, Poland\]
- [get_hf_thz.m](get_hf_thz.m) : Helper function to generate pseudo-random numbers from the FTR distribution. This is an efficient alternative to randpdf.m, particularly for FTR Monte-Carlo simulations. \[For details, refer to J. M. Romero-Jerez, F. J. Lopez-Martinez, J. F. Paris, and A. J. Goldsmith, “The fluctuating two-ray fading model: Statistical
characterization and performance analysis,” IEEE Trans. Wireless Commun., vol. 16, no. 7, pp. 4420–4432, Jul. 2017.\]
- [Combined_thz_experiment.m](Combined_thz_experiment.m) : A sample virtual experiment to demonstrate the use of all the above functions. Feel free to make changes as required for your own purpose.  

## Contributors:

Atharva Anand Joshi, Undergraduate Researcher 

Pranay Bhardwaj, PhD Researcher

Dr S. M. Zafaruddin, Senior Member, IEEE

Birla Institute of Technology and Science, Pilani


