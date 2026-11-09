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

(7) * PR Trenhago et al.
    * An integrated mathematical model of the cardiovascular and respiratory systems.
    * International Journal for Numerical Methods in Biomedical Engineering, 32(1):e02736, 2016. https://doi.org/10.1002/cnm.2736

(8) * TG Christensen and C Draeby.
    * Modelling the Respiratory System. MA thesis, Roskilde University, 1996.

## pH-extended version
This branch contains the pH-extended version of the model, based on version v2.0.0. The differences with respect to v2.0.0 are listed below.

Code:
- `transport.cc`, `transport.h`: the closed-form O2 and CO2 dissociation curves are replaced by a pH-dependent formulation following (7). Blood pH is computed by solving a cubic equation in [H+] with the Newton-Raphson method (`calculate_pH`), and CO2 is described by separate erythrocyte and plasma fractions. The curves are evaluated in mol/L and kPa and are inverted numerically (`invertedDissociation`). The initial conditions of gas transport are assigned as partial pressures.
- `closedloop.cc`: the arterial O2 saturation is obtained from the dissociation curves.
- `Makefile`: the executable is `cardiorespiratory_ph.x`.

Parameters: the parametrization is the same as in v2.0.0, except for the dissociation curve parameters.
- `transport.dat`: new section [pHmodel] with the parameters of the pH model and of the dissociation curves; [HPr]0 and [NaOH]0 are tuned, the other parameters are taken from (8). The initial conditions are partial pressures instead of concentrations, with the additional initial values H and pH.
- `commonParams_noControl.dat`: contains the steady-state values of the controlled variables of the pH-extended model.

Output:
- `GAS_AlveoliLungsAlg.txt`: [14] pHsysArt (systemic arterial pH), [15] pHpulCap (pulmonary capillary pH), [16] sao2 (arterial O2 saturation), [17] sppo2 (pulmonary capillary O2 saturation). In v2.0.0, column [14] is sao2.

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
- run "./cardiorespiratory_ph.x inputData.dat".

## Simulation settings (inputData.dat)

- [methods] method: numerical method, 0 for Euler, 1 for Heun, 2 for Runge-Kutta order 4.
- [convergence]: the simulation stops before "timeStop" when it has reached a steady state. From "nextConvCheck" on, every "windowConv" seconds, the mean values over the last window of systemic arterial pressure, central venous pressure, aortic valve flow, heart period (if the cardiovascular control is active), PaO2 and PaCO2 (if gas transport is active) are compared with the mean values over the previous window. If all relative changes are below "tolConv", the simulation stops and "convergence_report.txt" is written in the results folder. The check uses the sampled values, so sampling ("tSampleIni", "dtSample") must cover at least the last two windows. To disable the check, set "nextConvCheck" larger than "timeStop".
- [systems]: the systems are activated sequentially (cardiovascular system always active -> doLung -> doGas -> doControlCardio -> doControlResp). For example, it is not possible to activate doGas and deactivate doLung.
- [commonParams] commonParamsNoControl: if doControlCardio = 0, the common parameters are read from this file instead of "commonParams". "commonParams_noControl.dat" contains the controlled variables (resistances, unstressed volumes, maximum elastances, heart period) at steady state of a simulation with all controls active, so that a simulation without cardiovascular control reproduces the same baseline state.

If a NaN appears in the state variables, the simulation stops and "nan_detected.txt" is written in the results folder.
