/*
 * Copyright 2025 Bianca Maria Laudenzi, Caterina Dalmaso, Lucas Omar Muller
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

/**
 * @brief Gas transport and exchange model with pH-extended dissociation curves
 * 
 * @note Gas transport and exchange model proposed in:
 *
 * A Albanese, L Cheng, M Ursino & N. W. Chbat.
 * An integrated mathematical model of the human cardiopulmonary system: model development.
 * American Journal of Physiology-Heart and Circulatory Physiology 2016 310, H899-H921.   (A)
 *
 * @note O2 and CO2 dissociation curves and blood pH model following:
 *
 * P R Trenhago et al.
 * An integrated mathematical model of the cardiovascular and respiratory systems.
 * International Journal for Numerical Methods in Biomedical Engineering 2016 32(1), e02736.   (T)
 *
 * Units:
 * - concentration [mL/mL]
 * - pressure [mmHg]
 * - consumption rate [mL/min]
 * - volumes [mL]
 *
 * Units in the dissociation curves and pH model (calculate_pH, dissociation, invertedDissociation):
 * - concentration [mol/L] (from mL/mL: divide by 22.7 mL/mmol)
 * - pressure [kPa] (from mmHg: multiply by 0.133322 kPa/mmHg)
 * - hydrogen ion concentration [H+] [mol/L]
*/

class gasTransport
{
  public:
    /// time 
    double time;
    int verbose;
	int initial_cond = 1; // 0: assign concentrations, 1: assign pressures

	//-----------------------------------------------------------------------------------
    // State Variables
    //-----------------------------------------------------------------------------------


    /**
     * @brief Vector of 28 components that contains state variables
     *       fdo2	  [0]   // O2 dead space volume fraction 
     *       fdco2    [1]   // CO2 dead space volume fraction
     *       fAo2     [2]	// 02 alveoli volume fraction
     *       fAco2    [3]   // C02 alveoli volume fraction
	 *       cppo2	  [4]  // 02 pulmonary blood concentration 
	 *       cppco2   [5]  // C02 pulmonary blood concentration 
	 * // peripheral tissues gas concentrations 
	 *       cepo2	  [6]  // "e" extrasplanchic compartment
	 *       cepco2   [7]
	 *       cspo2	  [8]  // "s" splanchic compartment
	 *       cspco2   [9]
     *       cmpo2    [10]   // "m" skeletal muscle compartment
     *       cmpco2   [11]
	 *       chpo2    [12]   // "h" coronary compartment
	 *       chpco2   [13]
	 *       cbpo2    [14]   // "b" brain compartment
	 *       cbpco2   [15]
	 *  // venous pool gas concentrations
	 *       cevo2	  [16]   // "e" extrasplanchic compartment
	 *       cevco2   [17]
	 *       csvo2	  [18]   // "s" splanchic compartment
	 *       csvco2   [19]
	 *       cmvo2	  [20]   // "m" skeletal muscle compartment
	 *       cmvco2   [21]
     *       chvo2	  [22]   // "h" coronary compartment
     *       chvco2   [23]
     *       cbvo2	  [24]   // "b" brain compartment
     *       cbvco2   [25]
     *       cvo2	  [26]   // thoracic veins compartment
     *       cvco2    [27]
    */
    vector<double> stateVars = vector<double>(28);

	/**
	 * @brief Vector of 28 components that contains the derivate of the state variables with respect to time (dvdt)
	 *       dfdo2dt	  [0]   // O2 dead space volume fraction
	 *       dfdco2dt     [1]   // CO2 dead space volume fraction
	 *       dfAo2dt      [2]	// 02 alveoli volume fraction
	 *       dfAco2dt     [3]   // C02 alveoli volume fraction
	 *       dcppo2dt	  [4]  // 02 pulmonary blood concentration
	 *       dcppco2dt    [5]  // C02 pulmonary blood concentration 
	 * // peripheral tissues gas concentrations
	 *       dcepo2	      [6]  // "e" extrasplanchic compartment
	 *       dcepco2dt    [7]
	 *       dcspo2dt	  [8]  // "s" splanchic compartment
	 *       dcspco2dt    [9]
	 *       dcmpo2dt     [10]   // "m" skeletal muscle compartment
	 *       dcmpco2dt    [11]
	 *       dchpo2dt     [12]   // "h" coronary compartment
	 *       dchpco2dt    [13]
	 *       dcbpo2dt     [14]   // "b" brain compartment
	 *       dcbpco2dt    [15]
	 *  // venous pool gas concentrations
	 *       dcevo2dt	  [16]   // "e" extrasplanchic compartment
	 *       dcevco2dt    [17]
	 *       dcsvo2dt	  [18]   // "s" splanchic compartment
	 *       dcsvco2dt    [19]
	 *       dcmvo2dt	  [20]   // "m" skeletal muscle compartment
	 *       dcmvco2dt    [21]
	 *       dchvo2dt	  [22]   // "h" coronary compartment
	 *       dchvco2dt    [23]
	 *       dcbvo2dt	  [24]   // "b" brain compartment
	 *       dcbvco2dt    [25]
	 *       dcvo2dt	  [26]   // thoracic veins compartment
	 *       dcvco2dt     [27]
	*/
    vector<double> dvdt = vector<double>(28);
  
    //-----------------------------------------------------------------------------------
    // Parameters
    //-----------------------------------------------------------------------------------

    /**
	 * @brief Vector of 19 components that contains parameters for lung gas exchange and transport (tab. 4 - A)
     * // environmental conditions
	 *	fio2	[0]        // O2 fraction in inspired air
	 *	fico2	[1]        // CO2 farction in inspired air
	 *	k		[2]        // proportionality constant that allows convertion of volumes from body
	 *	patm	[3]        // atmospheric pressure
	 *	pws		[4]        // water vapour pressure
	 * // exchange coefficient that summarize the entire respiratory membrane properties (e.g. surface area, thickness, etc.) 
	 * kO2      [5]
	 * kCO2     [6]
    */
    vector<double> LungsGasParams = vector<double>(7);  

    /**
	 * @brief Vector of 15 components that contains parameters for tissue gas exchange and transport (Tab. 5 - A)
     * // "vt" tissue volume
     * // "mo2" metabolic oxygen consumption rate
     * // "mco2" metabolic co2 generation rate
	 *
	 * vtep	    [0]		// "e" extrasplanchic compartment
	 * mo2ep	[1]
	 * mco2ep	[2]
	 * vtsp	    [3]		// "s" splanchic compartment
	 * mo2sp	[4]
	 * mco2sp   [5]
	 * vtmp	    [6]		// "m" skeletal muscle compartment
	 * mo2mp	[7]
	 * mco2mp	[8]
     * vthp  	[9]  	// "h" coronary compartment 
     * mo2hp	[10]
     * mco2hp	[11]
     * vtbp	    [12]	// "b" brain compartment
     * mo2bp	[13]		
     * mco2bp	[14]
     */
    vector<double> TissuesVenousGasParam = vector<double>(15);   

	/**
	 * @brief Vector of 14 components that contains parameters for pH and dissociation curves
	 * [0] KaCO2 = pow(10, -6.1);      		// [mol/l] dissociation constant of CO2 (42)
	 * [1] KaPr = pow(10, -7.3);       		// [mol/l] dissociation constant of plasma proteins (43)
	 * [2] NaOH_0 = 43.0e-3;     		  	// [mol/l] [NaOH]  (46)
	 * [3] HPr_0 = 41.0e-3;      		  	// [mol/l] [HPr]  (47)
	 * [4] xHbf = 0.0; 						// [-] fraction of fetal hemoglobin concentration in the blood
	 * [5] xHbCO = 0.005; 					// [-] fraction of carboxyhemoglobin in the blood
	 * [6] xHi = 0.005; 					// [-] fraction of methemoglobin (hemiglobin) in the blood
	 * [7] cDPG = 5.01e-3 ; 				// [mol/l] concentration of 2,3-diphosphogycerate in the erythrocyte
	 * [8] T_C = 37.0; 						// [°C] blood temperature in Celsius
	 * [9] alphaO2 = 9.83e-6;  				// [mol/l/kPa] solubility of O2
	 * [10] cHb = 9.3e-3; 					// [mol/l] total concentration of Hb 
	 * [11] alphaCO2_Ery = 0.195e-3;		// [mol/l/kPA] solubility of CO2 in erythrocyte (33)
	 * [12] alphaCO2_Pla = 0.230e-3;		// [mol/l/kPA] solubility of CO2 in plasma (34)
	 * [13] cHb_Ery = 21.0e-3;				// [mol/l] concentration of Hb in erythrocyte 
	*/
	vector<double> GasParams = vector<double>(14);   

	/**
	 * @brief Vector of 2 components that contains parameters for 
	 * tault	    [0]		
	 * tauvl		[1]
     */
	vector<double> DelayedGasParam = vector<double>(2); 

	/**
	* @brief Vector of  components that contains parameters from other systems (cardiovascular, ) 
	* vd	[0]		// dead space volume       double vd
	* vdot	[1]		// total airflow           double Vdot
	* vadot [2]		// alveolar airflow        double VAdot
	* va	[3]		// volume in the alveoli   Lungs.vol[3]
	* vep	[4]		// CVS.stateVars[1-5]
	* vsp	[5]
	* vmp	[6]
	* vhp	[7]
	* vbp	[8]
	* vmv	[9]		// CVS.sysCircAlg[27]
	* vhv	[10]    // CVS.stateVars[9-11]
	* vbv	[11]
	* vtv	[12]    //                      + CVS.TveinParam[4]
	* vpp	[13]    // CVS.stateVars[13]
	* qep	[14]	// CVS.sysCircAlg[7-11]
	* qsp	[15]
	* qmp	[16]
	* qhp	[17]   
	* qbp	[18]
	* qev	[19]	// CVS.sysCircAlg[20-24]
	* qsv	[20]
	* qmv	[21]
	* qhv	[22]
	* qbv	[23] 
	* vev	[24]	// CVS.sysCircAlg[26-32]
	* vsv	[25]
	* qepin [26]
	* qspin [27]
	* qmpin [28]
	* qhpin [29]
	* qbpin [30]
	* qpp	[31]	//	pulCircAlg[2-3]
	* qps	[32]
	* qppin [33]
	*/
	vector<double> CommonVars = vector<double>(34);

    //-----------------------------------------------------------------------------------
    // Algebraic variables
    //-----------------------------------------------------------------------------------

    /**
	 * @brief Vector of 13 components containing the alveoli-lungs gas exchange
     * uvdot     [0]   // positive heaviside functions
	 * umvdot    [1]   // negative heaviside functions  
	 * pAco2     [2]   // CO2 partial pressure in the alveoli
     * pppo2     [3]   // O2 partial pressure in the pulmonary capillaries
	 * pppco2    [4]   // CO2 partial pressure in the pulmonary capillaries
	 * cao2      [5]   // O2 concentration in the systemic arteries
	 * caco2     [6]   // CO2 concentration in the systemic arteries 
	 * pao2      [7]   // O2 partial pressure in the systemic arteries
	 * paco2     [8]   // CO2 partial pressure in the systemic arteries
	 * cao2diss  [9]   // 
	 * mdoto2    [10]   // exchanged O2 between alveoli and capillaries
	 * mdotco2   [11]   // exchanged CO2 between alveoli and capillaries
	 * pAo2      [12]   // O2 partial pressure in the alveoli
	 */
	vector<double> AlveoliLungsAlg = vector<double>(13);

	/** 
	 * @brief Vector of 4 components containing the delayed systemic arterial/venous O2/CO2 concentrations
	 * cao2tilde  [0]
	 * caco2tilde [1]
	 * cvo2tilde  [2]
	 * cvco2tilde [3]
	 */
	vector<double> TissueVenousGasAlg = vector<double>(4);

	double sao2; // saturation systemic arteries
	double sppo2; // saturation pumonary capillary
	double HsysArt;
	double HpulCap;
	double pHsysArt;
	double pHpulCap;

    //-----------------------------------------------------------------------------------
    // Initial conditions of pressures
    //-----------------------------------------------------------------------------------
	// pppco2 = 38.9974;
	// pppo2 = 104.082;
	// pepco2 = 37.6065;
	// pepo2 = 67.5991;
	// pspco2 = 45.6982;
	// pspo2 = 30.6311;
	// pmpco2 = 40.7813;
	// pmpo2 = 38.0859;
	// phpco2 = 49.6231;
	// phpo2 = 40.5419;
	// pbpco2 = 43.9942;
	// pbpo2 = 32.7514;
	// pevco2 = 37.6065;
	// pevo2 = 67.5991;
	// psvco2 = 45.6982;
	// psvo2 = 30.6311;
	// pmvco2 = 40.7813;
	// pmvo2 = 38.0859;
	// phvco2 = 49.6231;
	// phvo2 = 40.5419;
	// pbvco2 = 43.9942;
	// pbvo2 = 32.7514;
	// pvco2 = 42.9239;
	// pvo2 = 39.0809;
	vector<double> PressuresInitialCond = vector<double>(24);


    //-----------------------------------------------------------------------------------
    // Outputs
    //-----------------------------------------------------------------------------------

    ofstream sampleGASstateVars;
    ofstream sampleGASalveoliLungsAlg;
	ofstream sampleGAStissueVenousGasAlg;
    string nameGAS;

    //-----------------------------------------------------------------------------------
    // Methods
    //-----------------------------------------------------------------------------------

    void init(string ifile, string ifileC, string outDir);
    void readGASparameters(string _file, string _fileC);
    void printGASparameters();
    void InitialCond(string _file);
	tuple<double, double, double> dissociation(double pCO2, double pO2, double cCO2_OLD, double H_old);
	tuple<double, double, double> invertedDissociation(double cCO2_target, double cO2_target, double H_old, double pCO2_OLD, double pO2_OLD);
	tuple<double, double> calculate_pH(double cCO2, double H_old);
	double calculate_SO2(double pO2, double pCO2, double pH);
	double calculate_cO2(double pO2, double SO2);
	double calculate_cCO2(double pCO2, double pH, double SO2);
	double df_dcCO2(double cCO2, double pCO2, double pO2, double H_old);
	double heaviside(double value1, double value2);
    void getAlgebraicRelations();
    vector<double> getTimeDerivative(vector<double> x);
	void additionalRelations();
	void updateFlowFilter(double dt);
	double qpp_filt = 0.;
	double qps_filt = 0.;
	double tauFlowFilter = 0.75;
	void getVariablesFromState(vector<double> x);
    void output();
};