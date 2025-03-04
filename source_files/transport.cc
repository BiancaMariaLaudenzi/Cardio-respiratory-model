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
	LungsGasParams[5] = ifile("LungsGas/csato2",0.);
	LungsGasParams[6] = ifile("LungsGas/csatco2",0.);
    LungsGasParams[7] = ifile("LungsGas/h1",0.);
    LungsGasParams[8] = ifile("LungsGas/h2",0.);
    LungsGasParams[9] = ifile("LungsGas/alpha1",0.);
    LungsGasParams[10] = ifile("LungsGas/alpha2",0.);
    LungsGasParams[11] = ifile("LungsGas/beta1",0.);
    LungsGasParams[12] = ifile("LungsGas/beta2",0.);
    LungsGasParams[13] = ifile("LungsGas/k1",0.);
    LungsGasParams[14] = ifile("LungsGas/k2",0.);
	LungsGasParams[15] = ifileC("shunt/sh", 0.);
	LungsGasParams[16] = ifile("LungsGas/hgb", 0.);
	LungsGasParams[17] = ifile("LungsGas/kO2", 0.);
	LungsGasParams[18] = ifile("LungsGas/kCO2", 0.);


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
	
	// [DelayedGasParam] (in closedloop.cc)
	DelayedGasParam[0] = ifile("DelayedGasParam/tault",0.);
	DelayedGasParam[1] = ifile("DelayedGasParam/tauvl",0.);

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


tuple<double, double> gasTransport::dissociation(double p1, double p2){
	/**
	* takes PO2(p1) and PCO2(p2) partial pressures and provides concentration of O2 and CO2 in blood,
	* as in Spencer(1979), eqs(4) and (5)
	*/
	double f1;
	double f2;
	double c1;
	double c2;

	f1 = p1 * (1. + LungsGasParams[11] * p2) / LungsGasParams[13] / (1. + LungsGasParams[9] * p2);
	f2 = p2 * (1. + LungsGasParams[12] * p1) / LungsGasParams[14] / (1. + LungsGasParams[10] * p1);

	c1 = LungsGasParams[5] * pow(f1,(1. / LungsGasParams[7])) / (1. + pow(f1,(1. / LungsGasParams[7])));
	c2 = LungsGasParams[6] * pow(f2,(1. / LungsGasParams[8])) / (1. + pow(f2,(1. / LungsGasParams[8])));

	return make_tuple(c1, c2);
}

tuple<double, double> gasTransport::invertedDissociation(double c1, double c2){
	/**
	* takes O2 and CO2 concentration and provides partial O2 and CO2 blood pressures,
	* as in Spencer(1979), eqs(6) and (7)
	*/
	double d1;
	double d2;
	double r1;
	double r2;
	double s1;
	double s2;
	double p1;
	double p2;

	d1 = max(LungsGasParams[13]*pow((c1 / (LungsGasParams[5] - c1)),LungsGasParams[7]), 0.);
	d2 = max(LungsGasParams[14]*pow((c2 / (LungsGasParams[6] - c2)),LungsGasParams[8]), 0.);

	r1 = -1.*(1. + LungsGasParams[11] * d2 - LungsGasParams[12] * d1 - LungsGasParams[9] * LungsGasParams[10]*d1*d2) / 2. / (LungsGasParams[12] + LungsGasParams[10] * LungsGasParams[11]*d2);
	r2 = -1.*(1. + LungsGasParams[12] * d1 - LungsGasParams[11] * d2 - LungsGasParams[10] * LungsGasParams[9]*d2*d1) / 2. / (LungsGasParams[11] + LungsGasParams[9] * LungsGasParams[12]*d1);

	s1 = -1.*(d1 + LungsGasParams[9] * d1*d2) / (LungsGasParams[12] + LungsGasParams[10] * LungsGasParams[11]*d2);
	s2 = -1.*(d2 + LungsGasParams[10] * d2*d1) / (LungsGasParams[11] + LungsGasParams[9] * LungsGasParams[12]*d1);

	p1 = r1 + pow((pow(r1,2) - s1),0.5);
	p2 = d2 * (1 + LungsGasParams[10] * p1) / (1 + LungsGasParams[12] * p1);
	
	return make_tuple(p1, p2);
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

	tie(AlveoliLungsAlg[3], AlveoliLungsAlg[4]) = invertedDissociation(stateVars[4], stateVars[5]); // pppo2,pppco2

	double denom = CommonVars[31] + CommonVars[32];
	if (denom != 0) {
		AlveoliLungsAlg[5] = max((CommonVars[31]*stateVars[4] + CommonVars[32]*TissueVenousGasAlg[2]) / denom, 0.); //A54 - cao2
		AlveoliLungsAlg[6] = max((CommonVars[31]*stateVars[5] + CommonVars[32]*TissueVenousGasAlg[3]) / denom, 0.); //A55 - caco2
	} else {
		AlveoliLungsAlg[5] = 0.;  
		AlveoliLungsAlg[6] = 0.;  
	}


	tie(AlveoliLungsAlg[7], AlveoliLungsAlg[8]) = invertedDissociation(AlveoliLungsAlg[5], AlveoliLungsAlg[6]); //pao2,paco2
	
	if ((LungsGasParams[5] - LungsGasParams[16]*1.34 - AlveoliLungsAlg[7]*0.003 / 100.) >= 0.) {
		AlveoliLungsAlg[9] = AlveoliLungsAlg[7] * 0.003 / 100.; //cao2diss
	}
	else {
		AlveoliLungsAlg[9] = LungsGasParams[5] - LungsGasParams[16] * 1.34;  //cao2diss
	}
	AlveoliLungsAlg[10] = LungsGasParams[17] * (AlveoliLungsAlg[12] - AlveoliLungsAlg[3]);  //mdoto2
	AlveoliLungsAlg[11] = LungsGasParams[18] * (AlveoliLungsAlg[2] - AlveoliLungsAlg[4]);  //mdotco2
}

void gasTransport::getTimeDerivative() {
	getAlgebraicRelations();
	// lung exchange
	dvdt[0] = 1. / CommonVars[0]*(AlveoliLungsAlg[0] * CommonVars[1]*(LungsGasParams[0] - max(stateVars[0], 0.)) + AlveoliLungsAlg[1] * CommonVars[2]*(max(stateVars[0], 0.) - max(stateVars[2], 0.))); // A42
	dvdt[1] = 1. / CommonVars[0]*(AlveoliLungsAlg[0] * CommonVars[1]*(LungsGasParams[1] - max(stateVars[1], 0.)) + AlveoliLungsAlg[1] * CommonVars[2]*(max(stateVars[1], 0.) - max(stateVars[3], 0.))); // A43
	dvdt[2] = 1. / CommonVars[3]*(AlveoliLungsAlg[0] * CommonVars[2]*(max(stateVars[0], 0.) - max(stateVars[2], 0.)) - LungsGasParams[2] * (AlveoliLungsAlg[10])); // A44
	dvdt[3] = 1. / CommonVars[3]*(AlveoliLungsAlg[0] * CommonVars[2]*(max(stateVars[1], 0.) - max(stateVars[3], 0.)) - LungsGasParams[2] * (AlveoliLungsAlg[11]));  // A45
	dvdt[4] = 1. / CommonVars[13]*(CommonVars[31]*(TissueVenousGasAlg[2] - stateVars[4]) + AlveoliLungsAlg[10]);
	dvdt[5] = 1. / CommonVars[13]*(CommonVars[31]*(TissueVenousGasAlg[3] - stateVars[5]) + AlveoliLungsAlg[11]);

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
}

void gasTransport::additionalRelations() {
	sao2 = (AlveoliLungsAlg[5] - AlveoliLungsAlg[9])/LungsGasParams[16]/1.34; // A56;
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
