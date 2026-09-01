# SEM-Optics
Defocus based phase retrieval of the electron probe wave function using focal series of probe intensities in SEM. 

#### Description #####

All functions and code are annotated and commented well and are clear to understand when read with the main paper. It takes a decent computational capability to run this program quickly.

main_code.m - Run this code to execute the program

ASM_propagation.m - Propagation function based on the angular spectrum method

defocus_phase_retrieval.m - Automated function to perform iterative phase retrieval based on defocus variation

volt2wavelen.m - Function to calculate the relativistic wavelength from accelerating

psf_visualization.m - Function to visualize the PSF of SEM Optics based on recovered probe wavefunction

z_sampling.mat - z-locations (mm) of the specimen plane in the focal series

probe_intensity_focal_series - Folder containing probe intensity images named (0_PSF.png, +-1.png, ...), although they are named as PSF.png, but they are not to be confused with PSF of the SEM Optics, they are just probe intensity.

**IMPORTANT DISTINCTION:** The probe intensity is loosely referred to as the PSF of SEM in SEM literature (inaccurate). We are visualizing the PSF of "SEM Optics - the lens system", which is a mathematically defined quantity.   

#### References: ####

**Main work:** Kamal, S., Nyaburi, K. M., & Hailstone, R. K. (2026). Optical point-spread function via defocus phase retrieval in scanning electron microscopy. Micron, 104097. https://doi.org/10.1016/j.micron.2026.104097

**Probe intensity reconstruction:** Zotta, M. D., Nevins, M. C., Hailstone, R. K., & Lifshin, E. (2018). The Determination and Application of the Point Spread Function in the Scanning Electron Microscope. Microscopy and microanalysis: the official journal of Microscopy Society of America, Microbeam Analysis Society, Microscopical Society of Canada, 24(4), 396–405. https://doi.org/10.1017/S1431927618012412

If you have any questions or suggestions about the program, please email: surya.kamal@rit.edu
