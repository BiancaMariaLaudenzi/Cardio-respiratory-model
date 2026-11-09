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

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include "GetPot.h"
#include <string>
#include <tuple>

#include "transport.h"

using namespace std;
using std::scientific;

void gasTransport::init(string ifile, string ifileC, string outDir) {

	// read parameter file
	if (verbose)
		cout << "Reading GAS params ..." << endl;
	if (ifile != "NOFILE") {
        if (ifileC !="NOFILE"){
		    readGASparameters(ifile,ifileC);
        }
	}
	if (verbose)
		cout << "Done reading GAS params" << endl;
		
	// intialiaze output file
	nameGAS = "GAS";
    string statenameStateVars = outDir + nameGAS + "_stateVars.txt";
	sampleGASstateVars.open(statenameStateVars.c_str());
    sampleGASstateVars << "# [0]: time;" << " "
					<< "[1]: fdo2;" << " "
                    << "[2]: fdco2;" << " "
                    << "[3]: fAo2;" << " "
                    << "[4]: fAco2;" << " "
                    << "[5]: cppo2;" << " "
                    << "[6]: cppco2;" << " "
                    << "[7]: cepo2;" << " "
                    << "[8]: cepco2;" << " "
                    << "[9]: cspo2;" << " "
                    << "[10]: cspco2;" << " "
                    << "[11]: cmpo2;" << " "
                    << "[12]: cmpco2;" << " "
                    << "[13]: chpo2;" << " " 
                    << "[14]: chpco2;" << " " 
                    << "[15]: cbpo2;" << " " 
                    << "[16]: cbpco2;" << " " 
					<< "[17]: cevo2;" << " " 
					<< "[18]: cevco2;" << " " 
                    << "[19]: csvo2;" << " " 
                    << "[20]: csvco2;" << " " 
                    << "[21]: cmvo2;" << " " 
				    << "[22]: cmvco2;" << " " 
                    << "[23]: chvo2;" << " " 
                    << "[24]: chvco2;" << " " 
                    << "[25]: cbvo2;" << " " 
                    << "[26]: cbvco2;" << " " 
                    << "[27]: cvo2;" << " " 
                    << "[28]: cvco2;" << " " 
                    << "\n";

	string statenameAlveoliLungsAlg = outDir + nameGAS + "_AlveoliLungsAlg.txt";
	sampleGASalveoliLungsAlg.open(statenameAlveoliLungsAlg.c_str());
	sampleGASalveoliLungsAlg << "# [0]: time;" << " "
                    << "[1]: uvdot;" << " "
                    << "[2]: umvdot;" << " "
                    << "[3]: pAco2;" << " "
                    << "[4]: pppo2;" << " "
                    << "[5]: pppco2;" << " "
                    << "[6]: cao2;" << " "
                    << "[7]: caco2;" << " "
                    << "[8]: pao2;" << " "
                    << "[9]: paco2;" << " "
                    << "[10]: cao2diss;" << " "
                    << "[11]: mdoto2;" << " "
                    << "[12]: mdotco2;" << " "
					<< "[13]: pAo2;" << " "
					<< "[14]: pHsysArt;" << " "
					<< "[15]: pHpulCap;" << " "
					<< "[16]: sao2;" << " "
					<< "[17]: sppo2;" << " "
					<< "\n";

	string statenameTissueVenousGasAlg = outDir + nameGAS + "_TissueVenousGasAlg(delayed).txt";
	sampleGAStissueVenousGasAlg.open(statenameTissueVenousGasAlg.c_str());
	sampleGAStissueVenousGasAlg << "# [0]: time;" << " "
                    << "[1]: cao2tilde;" << " "
                    << "[2]: caco2tilde;" << " "
                    << "[3]: cvo2tilde;" << " "
                    << "[4]: cvco2tilde;" << " "
					<< "\n";
					
	InitialCond(ifile);
	
	if (verbose){
		cout << "gasTransport::init :: DONE " << endl;
	}

}

void gasTransport::readGASparameters(string _file, string _fileC) {
	if (verbose)
		cout << " -- gasTransport::readGASparameters() reading from "<< _file << " and " << _fileC << endl;

	GetPot ifileC(_fileC.c_str());  // common parameters are assignes in commonParams.dat
	GetPot ifile(_file.c_str());
		
	// [Lungs Gas]
	LungsGasParams[0] = ifile("LungsGas/fio2",0.);
	LungsGasParams[1] = ifile("LungsGas/fico2",0.);
	LungsGasParams[2] = ifile("LungsGas/k",0.);
	LungsGasParams[3] = ifile("LungsGas/patm",0.);
	LungsGasParams[4] = ifile("LungsGas/pws", 0.);
	LungsGasParams[5] = ifile("LungsGas/kO2", 0.);
	LungsGasParams[6] = ifile("LungsGas/kCO2", 0.);

	// [Tissues Venous Gas] 
	// N.B. consumption rates have to be converted in mL/seconds.
	TissuesVenousGasParam[0] = ifile("TissuesVenousGas/vtep",0.);
	TissuesVenousGasParam[1] = ifile("TissuesVenousGas/mo2ep",0.)/60.;
    TissuesVenousGasParam[2] = ifile("TissuesVenousGas/mco2ep",0.)/60.;
	TissuesVenousGasParam[3] = ifile("TissuesVenousGas/vtsp", 0.);
    TissuesVenousGasParam[4] = ifile("TissuesVenousGas/mo2sp",0.)/60.;
    TissuesVenousGasParam[5] = ifile("TissuesVenousGas/mco2sp",0.)/60.;
    TissuesVenousGasParam[6] = ifile("TissuesVenousGas/vtmp",0.);
    TissuesVenousGasParam[7] = ifile("TissuesVenousGas/mo2mp",0.)/60.;
    TissuesVenousGasParam[8] = ifile("TissuesVenousGas/mco2mp",0.)/60.;
    TissuesVenousGasParam[9] = ifile("TissuesVenousGas/vthp",0.);
    TissuesVenousGasParam[10] = ifile("TissuesVenousGas/mo2hp",0.)/60.;
    TissuesVenousGasParam[11] = ifile("TissuesVenousGas/mco2hp",0.)/60.;
    TissuesVenousGasParam[12] = ifile("TissuesVenousGas/vtbp",0.);
    TissuesVenousGasParam[13] = ifile("TissuesVenousGas/mo2bp",0.)/60.;
    TissuesVenousGasParam[14] = ifile("TissuesVenousGas/mco2bp", 0.)/60.;
	
	// Gas Parameters
	GasParams[0] = ifile("pHmodel/KaCO2", pow(10., -6.1));      	
    GasParams[1] = ifile("pHmodel/KaPr", pow(10., -7.3));       	
    GasParams[2] = ifile("pHmodel/NaOH_0", 46.2e-3);     		
    GasParams[3] = ifile("pHmodel/HPr_0", 39.8e-3);      		
	GasParams[4] = ifile("pHmodel/xHbf", 0.0);                	
	GasParams[5] = ifile("pHmodel/xHbCO", 0.005);             	
	GasParams[6] = ifile("pHmodel/xHi", 0.005);               	
	GasParams[7] = ifile("pHmodel/cDPG", 5.01e-3);            
	GasParams[8] = ifile("pHmodel/T_C", 37.0);       	 		
	GasParams[9] = ifile("pHmodel/alphaO2", 9.83e-6);		
	GasParams[10] = ifile("pHmodel/cHb", 9.3e-3);			
	GasParams[11] = ifile("pHmodel/alphaCO2_Ery", 0.195e-3);	
	GasParams[12] = ifile("pHmodel/alphaCO2_Pla", 0.230e-3);	
	GasParams[13] = ifile("pHmodel/cHb_Ery", 21.0e-3);
	
	// [DelayedGasParam] (in closedloop.cc)
	DelayedGasParam[0] = ifile("DelayedGasParam/tault", 0.);
	DelayedGasParam[1] = ifile("DelayedGasParam/tauvl", 0.);

	if (verbose)
		printGASparameters();
}

void gasTransport::InitialCond(string _file) {
	if (verbose)
		cout << " -- gasTransport::saveInitialCond() reading from "
		<< _file << endl;

	GetPot ifile(_file.c_str());
	// [initialConditions]
	stateVars[0] = ifile("initial_conditions/fdo2", 0.);
	stateVars[1] = ifile("initial_conditions/fdco2", 0.);
	stateVars[2] = ifile("initial_conditions/fAo2", 0.);
	stateVars[3] = ifile("initial_conditions/fAco2", 0.);
	if (initial_cond == 0) {
		stateVars[4] = ifile("initial_conditions/cppo2", 0.);
		stateVars[5] = ifile("initial_conditions/cppco2", 0.);
		stateVars[6] = ifile("initial_conditions/cepo2", 0.);
		stateVars[7] = ifile("initial_conditions/cepco2", 0.);
		stateVars[8] = ifile("initial_conditions/cspo2", 0.);
		stateVars[9] = ifile("initial_conditions/cspco2", 0.);
		stateVars[10] = ifile("initial_conditions/cmpo2", 0.);
		stateVars[11] = ifile("initial_conditions/cmpco2", 0.);
		stateVars[12] = ifile("initial_conditions/chpo2", 0.);
		stateVars[13] = ifile("initial_conditions/chpco2", 0.);
		stateVars[14] = ifile("initial_conditions/cbpo2", 0.);
		stateVars[15] = ifile("initial_conditions/cbpco2", 0.);
		stateVars[16] = ifile("initial_conditions/cevo2", 0.);
		stateVars[17] = ifile("initial_conditions/cevco2", 0.);
		stateVars[18] = ifile("initial_conditions/csvo2", 0.);
		stateVars[19] = ifile("initial_conditions/csvco2", 0.);
		stateVars[20] = ifile("initial_conditions/cmvo2", 0.);
		stateVars[21] = ifile("initial_conditions/cmvco2", 0.);
		stateVars[22] = ifile("initial_conditions/chvo2", 0.);
		stateVars[23] = ifile("initial_conditions/chvco2", 0.);
		stateVars[24] = ifile("initial_conditions/cbvo2", 0.);
		stateVars[25] = ifile("initial_conditions/cbvco2", 0.);
		stateVars[26] = ifile("initial_conditions/cvo2", 0.);
		stateVars[27] = ifile("initial_conditions/cvco2", 0.);
	}
	HsysArt = ifile("initial_conditions/H", 0.);
	HpulCap = ifile("initial_conditions/H", 0.);
	pHsysArt = ifile("initial_conditions/pH", 0.);
	pHpulCap = ifile("initial_conditions/pH", 0.);
	AlveoliLungsAlg[7] = 90.0;
	AlveoliLungsAlg[8] = 40.0;
	double cCO2_old =  0.535371/22.7; // 22.0e-3; // mol/l

	if (initial_cond == 1) {
		PressuresInitialCond[0] = ifile("initial_conditions/pppco2", 38.9974);
		PressuresInitialCond[1] = ifile("initial_conditions/pppo2", 104.082);
		PressuresInitialCond[2] = ifile("initial_conditions/pepco2", 37.6065);
		PressuresInitialCond[3] = ifile("initial_conditions/pepo2", 67.5991);
		PressuresInitialCond[4] = ifile("initial_conditions/pspco2", 45.6982);
		PressuresInitialCond[5] = ifile("initial_conditions/pspo2", 30.6311);
		PressuresInitialCond[6] = ifile("initial_conditions/pmpco2", 40.7813);
		PressuresInitialCond[7] = ifile("initial_conditions/pmpo2", 38.0859);
		PressuresInitialCond[8] = ifile("initial_conditions/phpco2", 49.6231);
		PressuresInitialCond[9] = ifile("initial_conditions/phpo2", 40.5419);
		PressuresInitialCond[10] = ifile("initial_conditions/pbpco2", 43.9942);
		PressuresInitialCond[11] = ifile("initial_conditions/pbpo2", 32.7514);
		PressuresInitialCond[12] = ifile("initial_conditions/pevco2", 37.6065);
		PressuresInitialCond[13] = ifile("initial_conditions/pevo2", 67.5991);
		PressuresInitialCond[14] = ifile("initial_conditions/psvco2", 45.6982);
		PressuresInitialCond[15] = ifile("initial_conditions/psvo2", 30.6311);
		PressuresInitialCond[16] = ifile("initial_conditions/pmvco2", 40.7813);
		PressuresInitialCond[17] = ifile("initial_conditions/pmvo2", 38.0859);
		PressuresInitialCond[18] = ifile("initial_conditions/phvco2", 49.6231);
		PressuresInitialCond[19] = ifile("initial_conditions/phvo2", 40.5419);
		PressuresInitialCond[20] = ifile("initial_conditions/pbvco2", 43.9942);
		PressuresInitialCond[21] = ifile("initial_conditions/pbvo2", 32.7514);
		PressuresInitialCond[22] = ifile("initial_conditions/pvco2", 42.9239);
		PressuresInitialCond[23] = ifile("initial_conditions/pvo2", 39.0809);

		AlveoliLungsAlg[3] = PressuresInitialCond[1];
		AlveoliLungsAlg[4] = PressuresInitialCond[0];

		AlveoliLungsAlg[5] = 0.212045;  // mLgas/mLblood  
		AlveoliLungsAlg[6] = 0.535371;  // mLgas/mLblood

		for (int i = 0; i < 12; ++i) {
			double pCO2 = PressuresInitialCond[2 * i] * 0.133322;
			double pO2 = PressuresInitialCond[2 * i + 1] * 0.133322;
			double pH = -log10(HsysArt);
			double so2 = calculate_SO2(pO2, pCO2, HsysArt);
			tie(stateVars[2 * i + 4], stateVars[2 * i + 5], so2) = dissociation(pCO2, pO2, cCO2_old, HsysArt);
			cout << stateVars[2 * i + 4] << " " << stateVars[2 * i + 5] << endl;
			stateVars[2 * i + 4] *= 22.7;
			stateVars[2 * i + 5] *= 22.7;
		}
		}
	if (verbose){
		cout << "Printing initial conditions State Vars: ... " << endl;
		for (int i = 0; i < stateVars.size(); i++) {    
			cout << "State Vars" << i << " : " << stateVars[i] << endl;     
		}
	}
	
	getAlgebraicRelations();
}

void gasTransport::printGASparameters() {
	cout << "[Lung Gas Parama]" << endl;
	for (int i = 0; i < LungsGasParams.size(); i++) {
		cout << i << " value: " << LungsGasParams[i] << endl;
	}
	cout << "[Tissues Venous Gas Params]" << endl;
	for (int i = 0; i < TissuesVenousGasParam.size(); i++) {
		cout << i << " value: " << TissuesVenousGasParam[i] << endl;
	}
	cout << "delayed params: " << DelayedGasParam[0] << " " << DelayedGasParam[1] << endl;
}

tuple<double, double> gasTransport::calculate_pH(double cCO2, double H_old) {
	/**
	* Solves the cubic equation g([H+]) = 0 with Newton-Raphson and returns pH = -log10([H+]) 
    */
	const double rel_tol = 1e-10;     	  // tolerance
	const int max_iter = 100;

    // Coefficients of the cubic equation
    double a2 = GasParams[1] + GasParams[2] + GasParams[0];
    double a1 = GasParams[0] * (GasParams[2] - cCO2) + GasParams[1] * (GasParams[0] + GasParams[2] - GasParams[3]);
    double a0 = GasParams[0] * GasParams[1] * (GasParams[2] - GasParams[3] - cCO2);

    // Function g(H+) and its derivative
    auto g = [&](double H) {
        return pow(H, 3) + a2 * pow(H, 2) + a1 * H + a0;  // (48, 61)
    };
    auto dg = [&](double H) {
        return 3 * pow(H, 2) + 2 * a2 * H + a1;
    };

    // Newton–Raphson (62)
    double H = H_old;  // initial value [H+] mol/L
    for (int i = 0; i < max_iter; ++i) {
        double delta = g(H) / dg(H);
		double step = 1.0; // step size for damping
		double H_old_it = H;
		for (int j=0; j < max_iter ; ++j) {
			double H_new = H - step * delta;
			if (H_new > 0) { // ensure H remains positive
				H = H_new;
				break;
			}
			step *= 0.5;
		}
        if (fabs(H - H_old_it) < rel_tol * H_old_it) break;
    }
    return make_tuple(-log10(H), H);  // pH
}

double gasTransport::calculate_SO2(double pO2, double pCO2, double pH) {
	/**
     * Saturation of O2 in blood
     */
    double a = -0.72 * (pH - 7.4)										 // (27)
               + 0.09 * log(pCO2 / 5.33)
               + (0.07 - 0.03 * GasParams[4]) * (GasParams[7] - 5.0e-3)
               - 0.368 * GasParams[5]
               - 0.174 * GasParams[6]
               - 0.28 * GasParams[4];

    double x0 = 1.946 + a + 0.055 * (GasParams[8] - 37.0); 					 // (26)
    double x = log(pO2); 									 // (25)
    double h = 3.5 + a; 											 // (24)	
    double y = 1.875 + x - x0 + h * tanh(0.5343 * (x - x0));		 // (23)

    double SO2 = 1.0 / (1.0 + exp(-y)); 							 // (22)
	return SO2;
}


double gasTransport::calculate_cO2(double pO2, double SO2) {

	static bool first_call = true;
    /**
     * Total O2 concentration in blood = free + concentration in red blood cells
     */
	double cO2 = GasParams[9] * pO2 + GasParams[10] * SO2; 		// (21) 
    return cO2; 									// cO2 is returned in mol/l
}		 
	
double gasTransport::calculate_cCO2(double pCO2, double pH, double SO2) {
    /**
     *Concentration of CO2 in blood (erythrocytes + plasma)  
     */
	double cHbratio = GasParams[10]/GasParams[13]; 	// hemoglobin concentration of the blood (29)
	// pK erythrocyte and plasma
    double pHEry = 7.19 + 0.77 * (pH - 7.4) + 0.035 * (1.0 - SO2);					// (36)
	double pKEry = 6.125 - log10(1.0 + pow(10.0, pHEry - 7.84 - 0.06 * SO2)); 		// (35)
    double pKPla = 6.125 - log10(1.0 + pow(10.0, pH - 8.7)); //+ GasParams[12];		// (37)

    // cCO2 erythrocyte and plasma
    double cCO2_Ery = GasParams[11] * pCO2 * (1.0 + pow(10.0, pHEry - pKEry));		// (33)
    double cCO2_Pla = GasParams[12] * pCO2 * (1.0 + pow(10.0, pH - pKPla));			// (34)

    double cCO2 = (cCO2_Ery * cHbratio + cCO2_Pla * (1.0 - cHbratio)); 				// (29)
    return cCO2;						// cCO2 is returned in mol/l
}

double gasTransport::df_dcCO2(double cCO2, double pCO2, double pO2, double H_old) {
	auto phi = [&](double c) {
		double pH = std::get<0>(calculate_pH(c, H_old));
		double SO2 = calculate_SO2(pO2, pCO2, pH);
		return calculate_cCO2(pCO2, pH, SO2);
	};
    const double h = 1e-6 * cCO2;
    return (phi(cCO2 + h) - phi(cCO2 - h)) / (2 * h) - 1.0;
}

tuple<double, double, double> gasTransport::dissociation(double pCO2, double pO2, double cCO2_old, double H_old) { 
    const double tol = 1e-4 * cCO2_old;
    const int max_iter = 100;

    double cCO2 = cCO2_old; // 22.0e-3 mol/l = 0.4994 mLgas/mLblood;  // initial

    for (int i = 0; i < max_iter; ++i) {
        double pH = get<0>(calculate_pH(cCO2, H_old));
		double SO2 = calculate_SO2(pO2, pCO2, pH);
        double phi = calculate_cCO2(pCO2, pH, SO2);
	double f = phi - cCO2;
	double df = df_dcCO2(cCO2, pCO2, pO2, H_old);
		
		// Newton-Raphson step
		double delta = f / df;
        cCO2 -= delta;

        if (fabs(delta) < tol) break;
    }

    double pH_out = get<0>(calculate_pH(cCO2, H_old));
	double SO2_out = calculate_SO2(pO2, pCO2, pH_out);
	double cO2_out = calculate_cO2(pO2, SO2_out); 
    double cCO2_out = cCO2;

	return make_tuple(cO2_out, cCO2_out, SO2_out);
}

tuple<double, double, double> gasTransport::invertedDissociation(double cO2_target_mL, double cCO2_target_mL, double H_old, double pO2_OLD, double pCO2_OLD) {

    double pCO2_mmHg = pCO2_OLD;   // initial guess
    double pO2_mmHg = pO2_OLD;     // initial guess

	// convert concentrations from mLgas/mLblood to mol/l
	double cCO2_target = cCO2_target_mL / 22.7; 
	double cO2_target = cO2_target_mL / 22.7;
	double cCO2_old = cCO2_target; // initial guess
	double cO2_old = cO2_target; // initial guess
	// convert pressures from mmHg to kPa 
	double pO2 =  pO2_mmHg * 0.133322; 
	double pCO2 =  pCO2_mmHg * 0.133322;

	// const double hCO2 = 1e-4 * pCO2;
	// const double hO2 = 1e-4 * pO2;
	const double hCO2 = std::max(1e-6, 1e-4 * fabs(pCO2));
	const double hO2  = std::max(1e-6, 1e-4 * fabs(pO2));
    const double tol = 1e-6;
    const int max_iter = 100;
	double cO2, cCO2;
	double cO2_pluspO2, cCO2_pluspO2;
	double cO2_minuspO2, cCO2_minuspO2;
	double cO2_pluspCO2, cCO2_pluspCO2;
	double cO2_minuspCO2, cCO2_minuspCO2;
	double SO2, SO2_pluspO2, SO2_minuspO2, SO2_pluspCO2, SO2_minuspCO2;

    for (int iter = 0; iter < max_iter; ++iter) {
    	tie(cO2, cCO2, SO2) = dissociation(pCO2, pO2, cCO2_old, H_old);
		tie(cO2_pluspO2, cCO2_pluspO2, SO2_pluspO2) = dissociation(pCO2, pO2 + hO2, cCO2_old, H_old);
		tie(cO2_minuspO2, cCO2_minuspO2, SO2_minuspO2) = dissociation(pCO2, pO2 - hO2, cCO2_old, H_old);
		tie(cO2_pluspCO2, cCO2_pluspCO2, SO2_pluspCO2) = dissociation(pCO2 + hCO2, pO2, cCO2_old, H_old);
		tie(cO2_minuspCO2, cCO2_minuspCO2, SO2_minuspCO2) = dissociation(pCO2 - hCO2, pO2, cCO2_old, H_old);

		double f1 = cO2 - cO2_target;
        double f2 = cCO2 - cCO2_target;

        // Numerical partial derivatives
        double df1_pO2   = (cO2_pluspO2 - cO2_minuspO2) / (2*hO2);
        double df1_pCO2  = (cO2_pluspCO2 - cO2_minuspCO2) / (2*hCO2);
        double df2_pO2  = (cCO2_pluspO2 - cCO2_minuspO2) / (2*hO2);
        double df2_pCO2 = (cCO2_pluspCO2 - cCO2_minuspCO2) / (2*hCO2);

        // Jacobian matrix
        double det = df2_pCO2 * df1_pO2 - df2_pO2 * df1_pCO2;
		if (fabs(det) < 1e-12) {
			cerr << "Jacobian near singular!" << endl;
			break;
		}

        // Inverse of the Jacobian * [f1, f2]
        double dpO2  = (f1 * df2_pCO2 / det - f2 * df1_pCO2 / det);
        double dpCO2 = (-f1 * df2_pO2 / det + f2 * df1_pO2 / det);

		// pCO2 -= dpCO2;
		// pO2  -= dpO2;

		double alpha = 0.5;
		pCO2 -= alpha * dpCO2;
		pO2  -= alpha * dpO2;

		double rel_err = std::max(
			abs(f1) / std::max(1e-8, fabs(cO2_target)),
			abs(f2) / std::max(1e-8, fabs(cCO2_target))
		);

        // double rel_err = std::max(abs(f1) / cO2_target, abs(f2) / cCO2_target);
		if (rel_err < tol)  break;
		if (iter == max_iter - 1) {
			cerr << "Warning: Inverted dissociation did not converge." << endl;}
    }

	// convert pressures from kPa to mmHg
	double pO2_out =  pO2 / 0.133322; 
	double pCO2_out =  pCO2 / 0.133322;

	return make_tuple(pCO2_out, pO2_out, SO2);
}

double gasTransport::heaviside(double value1, double value2){
	double value;
	if (value1 < 0) {
		value = 0;
	}
	if (value1 == 0) {
		value = value2;
	}
	if (value1 > 0) {
		value = 1;
	}
	return value;
}

void gasTransport::getAlgebraicRelations() {
	AlveoliLungsAlg[0] = heaviside(CommonVars[1], 0.);      //heaviside from vdot
	AlveoliLungsAlg[1] = heaviside(-1.*CommonVars[1], 0.);
		
	AlveoliLungsAlg[12] = max(stateVars[2], 0.)*(LungsGasParams[3] - LungsGasParams[4]);  // A - pAo2
	AlveoliLungsAlg[2] = max(stateVars[3], 0.)*(LungsGasParams[3] - LungsGasParams[4]); // A53 - pAco2

	if (stateVars[5] <= 0.01) {
		cout << "cppCO2 negative, setting to zero" << endl;
		stateVars[5] = 0.0001;
	}
	tie (pHpulCap, HpulCap) = calculate_pH(stateVars[5] / 22.7, HpulCap); // pH	
	tie(AlveoliLungsAlg[4], AlveoliLungsAlg[3], sppo2) = invertedDissociation(stateVars[4], stateVars[5], HpulCap, AlveoliLungsAlg[3], AlveoliLungsAlg[4]); // pppo2,pppco2

    double denom_filt = qpp_filt + qps_filt;
    if (denom_filt > 0.001) {
        AlveoliLungsAlg[5] = (qpp_filt * stateVars[4] + qps_filt * TissueVenousGasAlg[2]) / denom_filt; 
        AlveoliLungsAlg[6] = (qpp_filt * stateVars[5] + qps_filt * TissueVenousGasAlg[3]) / denom_filt; 
    }
	tie (pHsysArt, HsysArt) = calculate_pH(AlveoliLungsAlg[6]/ 22.7, HsysArt);
	tie(AlveoliLungsAlg[8], AlveoliLungsAlg[7], sao2) = invertedDissociation(AlveoliLungsAlg[5], AlveoliLungsAlg[6], HsysArt, AlveoliLungsAlg[7], AlveoliLungsAlg[8]); //pao2,paco2
	AlveoliLungsAlg[10] = LungsGasParams[5] * (AlveoliLungsAlg[12] - AlveoliLungsAlg[3]);  //mdoto2
	AlveoliLungsAlg[11] = LungsGasParams[6] * (AlveoliLungsAlg[2] - AlveoliLungsAlg[4]);  //mdotco2
}

vector<double> gasTransport::getTimeDerivative(vector<double> x) {
	for (int i = 0; i < stateVars.size(); i++) {
		stateVars[i] = x[i];
	}
	AlveoliLungsAlg[0] = heaviside(CommonVars[1], 0.);
	AlveoliLungsAlg[1] = heaviside(-1.*CommonVars[1], 0.);

	// lung exchange
	dvdt[0] = 1. / CommonVars[0]*(AlveoliLungsAlg[0] * CommonVars[1]*(LungsGasParams[0] - max(stateVars[0], 0.)) + AlveoliLungsAlg[1] * CommonVars[2]*(max(stateVars[0], 0.) - max(stateVars[2], 0.))); // A42
	dvdt[1] = 1. / CommonVars[0]*(AlveoliLungsAlg[0] * CommonVars[1]*(LungsGasParams[1] - max(stateVars[1], 0.)) + AlveoliLungsAlg[1] * CommonVars[2]*(max(stateVars[1], 0.) - max(stateVars[3], 0.))); // A43
	dvdt[2] = 1. / CommonVars[3]*(AlveoliLungsAlg[0] * CommonVars[2]*(max(stateVars[0], 0.) - max(stateVars[2], 0.)) - LungsGasParams[2] * (AlveoliLungsAlg[10])); // A44
	dvdt[3] = 1. / CommonVars[3]*(AlveoliLungsAlg[0] * CommonVars[2]*(max(stateVars[1], 0.) - max(stateVars[3], 0.)) - LungsGasParams[2] * (AlveoliLungsAlg[11]));  // A45
	dvdt[4] = 1. / CommonVars[13]*(CommonVars[33]*(TissueVenousGasAlg[2] - stateVars[4]) + AlveoliLungsAlg[10]);
	dvdt[5] = 1. / CommonVars[13]*(CommonVars[33]*(TissueVenousGasAlg[3] - stateVars[5]) + AlveoliLungsAlg[11]);

	// tissue gas exchange
	dvdt[6] = 1. / (TissuesVenousGasParam[0] + CommonVars[4])*(CommonVars[26]*(TissueVenousGasAlg[0] - max(stateVars[6], 0.)) - TissuesVenousGasParam[1]); // A63
	dvdt[7] = 1. / (TissuesVenousGasParam[0] + CommonVars[4])*(CommonVars[26]*(TissueVenousGasAlg[1] - max(stateVars[7], 0.)) + TissuesVenousGasParam[2]); // A64
	dvdt[8] = 1. / (TissuesVenousGasParam[3] + CommonVars[5])*(CommonVars[27]*(TissueVenousGasAlg[0] - max(stateVars[8], 0.)) - TissuesVenousGasParam[4]); // A65
	dvdt[9] = 1. / (TissuesVenousGasParam[3] + CommonVars[5])*(CommonVars[27]*(TissueVenousGasAlg[1] - max(stateVars[9], 0.)) + TissuesVenousGasParam[5]);  // A66
	dvdt[10] = 1. / (TissuesVenousGasParam[6] + CommonVars[6])*(CommonVars[28]*(TissueVenousGasAlg[0] - max(stateVars[10], 0.)) - TissuesVenousGasParam[7]); // A61
	dvdt[11] = 1. / (TissuesVenousGasParam[6] + CommonVars[6])*(CommonVars[28]*(TissueVenousGasAlg[1] - max(stateVars[11], 0.)) + TissuesVenousGasParam[8]); // A62
	dvdt[12] = 1. / (TissuesVenousGasParam[9] + CommonVars[7])*(CommonVars[29]*(TissueVenousGasAlg[0] - max(stateVars[12], 0.)) - TissuesVenousGasParam[10]); // A57
	dvdt[13] = 1. / (TissuesVenousGasParam[9] + CommonVars[7])*(CommonVars[29]*(TissueVenousGasAlg[1] - max(stateVars[13], 0.)) + TissuesVenousGasParam[11]); // A58
	dvdt[14] = 1. / (TissuesVenousGasParam[12] + CommonVars[8])*(CommonVars[30]*(TissueVenousGasAlg[0] - max(stateVars[14], 0.)) - TissuesVenousGasParam[13]); // A59
	dvdt[15] = 1. / (TissuesVenousGasParam[12] + CommonVars[8])*(CommonVars[30]*(TissueVenousGasAlg[1] - max(stateVars[15], 0.)) + TissuesVenousGasParam[14]); // A60

	// venous pool gas transport
	dvdt[16] = 1. / CommonVars[24]*CommonVars[14]*(max(stateVars[6], 0.) - max(stateVars[16], 0.)); // A73
	dvdt[17] = 1. / CommonVars[24]*CommonVars[14]*(max(stateVars[7], 0.) - max(stateVars[17], 0.)); // A74 
	dvdt[18] = 1. / CommonVars[25]*CommonVars[15]*(max(stateVars[8], 0.) - max(stateVars[18], 0.)); // A75
	dvdt[19] = 1. / CommonVars[25]*CommonVars[15]*(max(stateVars[9], 0.) - max(stateVars[19], 0.)); // A76
	dvdt[20] = 1. / CommonVars[9]*CommonVars[16]*(max(stateVars[10], 0.) - max(stateVars[20], 0.)); // A71
	dvdt[21] = 1. / CommonVars[9]*CommonVars[16]*(max(stateVars[11], 0.) - max(stateVars[21], 0.)); // A72 
	dvdt[22] = 1. / CommonVars[10]*CommonVars[17]*(max(stateVars[12], 0.) - max(stateVars[22], 0.)); // A67
	dvdt[23] = 1. / CommonVars[10]*CommonVars[17]*(max(stateVars[13], 0.) - max(stateVars[23], 0.)); // A68
	dvdt[24] = 1. / CommonVars[11]*CommonVars[18]*(max(stateVars[14], 0.) - max(stateVars[24], 0.));// A69
	dvdt[25] = 1. / CommonVars[11]*CommonVars[18]*(max(stateVars[15], 0.) - max(stateVars[25], 0.));// A70 

	dvdt[26] = 1. / (CommonVars[12])* (CommonVars[22]*(max(stateVars[22], 0.) - max(0., stateVars[26])) + CommonVars[23]*(max(stateVars[24], 0.) - max(0., stateVars[26])) + CommonVars[21]*(max(stateVars[20], 0.) - max(0., stateVars[26])) + CommonVars[19]*(max(stateVars[16], 0.) - max(0., stateVars[26])) + CommonVars[20]*(max(stateVars[18], 0.) - max(0., stateVars[26]))); // A77
	dvdt[27] = 1. / (CommonVars[12])* (CommonVars[22]*(max(stateVars[23], 0.) - max(0., stateVars[27])) + CommonVars[23]*(max(stateVars[25], 0.) - max(0., stateVars[27])) + CommonVars[21]*(max(stateVars[21], 0.) - max(0., stateVars[27])) + CommonVars[19]*(max(stateVars[17], 0.) - max(0., stateVars[27])) + CommonVars[20]*(max(stateVars[19], 0.) - max(0., stateVars[27])));  // A78

	return dvdt;
}

void gasTransport::getVariablesFromState(vector<double> x) {
    for (int i = 0; i < stateVars.size(); i++) {
        stateVars[i] = x[i];
    }
    getAlgebraicRelations();
}

void gasTransport::output(){
	sampleGASstateVars << scientific << setprecision(12);
    sampleGASstateVars << time << " ";
	// State variables
	for (int i = 0; i < stateVars.size(); i++) {    
		sampleGASstateVars << stateVars[i] << " ";     
	}
	sampleGASstateVars << "\n";
	sampleGASstateVars.flush();
	// Algebraic variables
	sampleGASalveoliLungsAlg << scientific << setprecision(12);
	sampleGASalveoliLungsAlg << time << " ";
	for (int i = 0; i < AlveoliLungsAlg.size(); i++) {
		sampleGASalveoliLungsAlg << AlveoliLungsAlg[i] << " ";
	}
	sampleGASalveoliLungsAlg << pHsysArt << " ";
	sampleGASalveoliLungsAlg << pHpulCap << " ";
	sampleGASalveoliLungsAlg << sao2 << " ";
	sampleGASalveoliLungsAlg << sppo2 << " ";
	sampleGASalveoliLungsAlg << "\n";
	sampleGASalveoliLungsAlg.flush();
	
	sampleGAStissueVenousGasAlg << scientific << setprecision(12);
	sampleGAStissueVenousGasAlg << time << " ";
	for (int i = 0; i < TissueVenousGasAlg.size(); i++) {
		sampleGAStissueVenousGasAlg << TissueVenousGasAlg[i] << " ";
	}
	sampleGAStissueVenousGasAlg << "\n";
	sampleGAStissueVenousGasAlg.flush();
}

void gasTransport::updateFlowFilter(double dt) {
	double qpp = fmax(CommonVars[31],0.);
	double qps = fmax(CommonVars[32],0.);
	double alpha = dt / (dt + tauFlowFilter);
	if (qpp_filt == 0.0 && qpp != 0.0) qpp_filt = qpp;
	if (qps_filt == 0.0 && qps != 0.0) qps_filt = qps;
	qpp_filt = alpha * qpp + (1.0 - alpha) * qpp_filt;
	qps_filt = alpha * qps + (1.0 - alpha) * qps_filt;
}
