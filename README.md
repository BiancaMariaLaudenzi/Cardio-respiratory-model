## License
This project is licensed under the **Apache License 2.0**. 
You can find the full text of the license in the `LICENSE` file or at the following link: 
🔗 https://www.apache.org/licenses/LICENSE-2.0 

## Third-Party Components
The file `GetPot.h` (modified from the original version, which can be found at https://sourceforge.net/projects/getpot/) is distributed under the **GNU LGPL v2.1** license. 
For more details, see `LICENSE.LGPL` or visit: 
🔗 https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html

## Data Sources
The data in the directory "params" is partially sourced from published scientific papers. The data is not covered by the Apache License 2.0 but remain under the copyright or license specified by their respective authors.

Sources:

(0) * BM Laudenzi, A Cucino, S Lassola, E Balzani, LO Muller
    * Predicting in-hospital indicators from wearable-derived signals for cardiovascular and respiratory disease monitoring: an in silico study.
    * Submitted for publication, March 2025

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

## How to run the code

The directory "source_files" contains the files ".cc" and ".h" representing the different systems: the cardiovascular system (cardiovascular), the lungs mechanics system (lungMechanics), the gas exchange and transport system (transport), the cardiovascular control system (cardiovascularControl) and the ventilatory control system (respiratoryControl). 

The "test/params" directory contains the input parameters of the model divided for each system. Note that the file "commonParams.pot" contains the input parameters that are common to more than one system. Results will be saved in a subfolder of the "test" directory. 

To run the code:
- enter the "source_files" directory and run "make -f Makefile clean all";
- exit the "source_files" diretory, enter the "test" directory and run "make -f Makefile clean all";
- define appropriate input data in "inputData.dat";
- run "./cardiorespiratory.x inputData.dat".


