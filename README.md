## License
This project is licensed under the **Apache License 2.0**. 
You can find the full text of the license in the `LICENSE` file or at the following link: 
🔗 https://www.apache.org/licenses/LICENSE-2.0 

## Third-Party Components
The file `GetPot.h` (modified from the original version, which can be found at https://sourceforge.net/projects/getpot/) is distributed under the **GNU LGPL v2.1** license. 
For more details, see `LICENSE.LGPL` or visit: 
🔗 https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html

## Data Sources
The data in the directory "baseline" is partially sourced from published scientific papers. The data is not covered by the Apache License 2.0 but remain under the copyright or license specified by their respective authors.

Sources:

(0) * BM Laudenzi, A Cucino, S Lassola, E Balzani, LO Müller
    * Predicting in-hospital indicators from wearable-derived signals for cardiovascular and respiratory disease monitoring: An in silico study.
    * PLOS Digital Health 4(10): e0001041, 2025. https://doi.org/10.1371/journal.pdig.0001041

(1) * A Albanese, L Cheng, M Ursino & N. W. Chbat.
    * An integrated mathematical model of the human cardiopulmonary system: model development.
    * American Journal of Physiology-Heart and Circulatory Physiology 2016 310, H899-H921. 

(2) * JP Mynard, MR Davidson, DJ Penny, and JJ Smolich. 
    * A simple, versatile valve model for use in the lumped parameter and one-dimensional cardiovascular models. 
    * International Journal for Numerical Methods in the Biomedical Engineering, 28(6-7):626–641, 2012. 

(3) * Mauro Ursino and Elisa Magosso. 
    * Acute cardiovascular response to isocapnic hypoxia. I. A mathematical model. 
    * American Journal of Physiology-Heart and Circulatory Physiology, 279(1):H149-H165, July 2000. 

(4) * FY Liang, S Takagi, R Himeno, and H Liu. 
    * Biomechanical characterization of ventricular–arterial coupling during aging: a multi-scale model study. 
    * Journal of biomechanics, 42(6):692–704, 2009. 

(5) * E. Magosso and M. Ursino. A mathematical model of CO2 effect on cardiovascular regulation.
    *   American Journal of Physiology. Heart and Circulatory Physiology, 281(5):H2036–2052, November 2001.
 
(6) * Mauro Ursino and Elisa Magosso. A theoretical analysis of the carotid body chemoreceptor
    *    response to O2 and CO2 pressure changes. Respiratory Physiology & Neurobiology, 130(1):99–110, March 2002.

## Changes
The changes with respect to the previous version are listed in `CHANGES.txt`.

## Directory structure

- "source_files": files ".cc" and ".h" representing the different systems: the cardiovascular system (cardiovascular), the lungs mechanics system (lungMechanics), the gas exchange and transport system (transport), the cardiovascular control system (cardiovascularControl) and the ventilatory control system (respiratoryControl).
- "baseline": input parameters of the model divided for each system. The file "commonParams.dat" contains the input parameters that are common to more than one system. The file "commonParams_noControl.dat" contains the common parameters used when the cardiovascular control is not active (see below).
- "inputData.dat": simulation settings (numerical method, time parameters, active systems, results folder, parameter files).
- "main.cc" and "Makefile": main program and Makefile to build the executable.

Results are saved in the folder defined by "pathResults" in "inputData.dat".

Please check that parameter files only use "LF" end-of-line characters. If end-of-line characters are "CRLF", parameters are not read correctly.

## How to run the code

- enter the "source_files" directory and run "make -f Makefile clean all";
- go back to the main directory and run "make -f Makefile clean all";
- define appropriate input data in "inputData.dat";
- run "./cardiorespiratory.x inputData.dat".

## Simulation settings (inputData.dat)

- [methods] method: numerical method, 0 for Euler, 1 for Heun, 2 for Runge-Kutta order 4.
- [convergence]: the simulation stops before "timeStop" when it has reached a steady state. From "nextConvCheck" on, every "windowConv" seconds, the mean values over the last window of systemic arterial pressure, central venous pressure, aortic valve flow, heart period (if the cardiovascular control is active), PaO2 and PaCO2 (if gas transport is active) are compared with the mean values over the previous window. If all relative changes are below "tolConv", the simulation stops and "convergence_report.txt" is written in the results folder. The check uses the sampled values, so sampling ("tSampleIni", "dtSample") must cover at least the last two windows. To disable the check, set "nextConvCheck" larger than "timeStop".
- [systems]: the systems are activated sequentially (cardiovascular system always active -> doLung -> doGas -> doControlCardio -> doControlResp). For example, it is not possible to activate doGas and deactivate doLung.
- [commonParams] commonParamsNoControl: if doControlCardio = 0, the common parameters are read from this file instead of "commonParams". "commonParams_noControl.dat" contains the controlled variables (resistances, unstressed volumes, maximum elastances, heart period) at steady state of a simulation with all controls active, so that a simulation without cardiovascular control reproduces the same baseline state.

If a NaN appears in the state variables, the simulation stops and "nan_detected.txt" is written in the results folder.
