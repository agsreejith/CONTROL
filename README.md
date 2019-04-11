# CONTROL
 This software is intented to be fully automated, aimed at producing science-quality output with a single command line with zero user interference for CUTE data. It can be easily used for any single order spectral data in any wavelength without any modification. 
 
The program works in a series of steps following standard CCD reduction techniques user can also select individual modules using the step function described above.
A reduction log is created to help users with processes carried out and mentioning the different parameters used if any.
The steps involved in the program execution are as follows:
 1. Prepare: Check for the different input parameters are set variables accordingly. If data files are present classify them according to     file type.
 2. Hot and bad pixel correction: Correct for hot and bad pixels in the frames.
 3. Create master bias: Create a master bias if master bias file is not found from a set of bias files.
	4. Create master dark: Create a master dark if master dark file is not found from a set of dark files.
	5. Create master flat: Create a master flat if master flat file is not found from a set of flat files.
	6. Correct for bias: Correct of the effect of bias in science frames.
	7. Correct for dark: Correct of the effect of dark in science frames.
	8. Correct for flat: Correct of the effect of flat in science frames.
	9. Correct for cosmic rays: Correct for the effect of cosmic rays using LA cosmic algorithm. 
	10.Extract spectrum: Define spectrum trace and extract the spectrum.
	11.Subtract Background: Subtract background from spectrum.
	12.Wavelength calibration: Do wavelength calibration
	13.Flux calibration: Do flux calibration
	14.Default Light curve: Create 3 light curves (short, middle and long wavelength)
	15.Transmission spectrum: Retrie transmission spectrum.
