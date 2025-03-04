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
#include <iostream>
#include <sstream>
#include <string>

#include "lungMechanics.h"

using namespace std;
using std::scientific;

void lungMechanics::init(string ifile, string ifileC, string outDir) {

	// read parameter file
	if (verbose)
		cout << "Reading lung mechanics params ..." << endl;
	if (ifile != "NOFILE") {
		if (ifileC !="NOFILE"){
			readLungMechParameters(ifile, ifileC);
		}
	}
	if (verbose)
		cout << "Done reading lung mechanics params" << endl;
	
	// intialiaze output file
	nameLung = "lungMech";
	string statenameVol = outDir + nameLung + "_Vol.txt";
	sampleLungVol.open(statenameVol.c_str());
	sampleLungVol << "# [0] : time:" << " "
					<< "[1] : vl;" << " " 
					<< "[2] : vtr;" << " "
					<< "[3] : vb;" << " "
					<< "[4] : vA;" << " "
					<< "[5] : vpl;" << " "
					<< "[6] : xi;" << " "
					<< "[7] : vd;" << " "
					<< "[8] : vtot" << " "
					<< "\n";

	string statenamePr = outDir + nameLung + "_Pr.txt";
	sampleLungPr.open(statenamePr.c_str());
	sampleLungPr << "# [0]: time;" << " "
					<< "[1]: pl;" << " " 
					<< "[2]: ptr;" << " "
					<< "[3]: pb;" << " "
					<< "[4]: pA;" << " "
					<< "[5]: ppl1;" << " "
					<< "[6]: Pmus"
					<< "[7]: pabd"
					<< "\n";

	string statenameFlow = outDir + nameLung + "_Flows.txt";
	sampleLungFlow.open(statenameFlow.c_str());
	sampleLungFlow << "# [0]: time;" << " "
					<< "[1]: Vdot;" << " " 
					<< "[2]: VAdot;"
					<< "\n";
	
	InitialCond();

	if (verbose){
		cout << "lungMechanics::init :: DONE " << endl;
	}

}


void lungMechanics::readLungMechParameters(string _file, string _fileC) {
	if (verbose)
		cout << " -- lungMechanics::readLungMechParameters() reading from "
				<< _file << " and "
                << _fileC << endl;

	GetPot ifile(_file.c_str());
	GetPot ifileC(_fileC.c_str());
	
	// [compliances]
	C[0] = ifile("compliances/cl", 0.00127); 
	C[1] = ifile("compliances/ctr", 0.00238);
	C[2] = ifile("compliances/cb", 0.0131);
	C[3] = ifile("compliances/cA", 0.2);
	C[4] = ifile("compliances/ccw", 0.2445);
	C[5] = ifile("compliances/cabd", 0.183);

	// [unstressed volumes]
	Vu[0] = ifile("unstressed_volumes/vul", 34.4);
	Vu[1] = ifile("unstressed_volumes/vutr", 6.63);
	Vu[2] = ifile("unstressed_volumes/vub", 18.7);
	Vu[3] = ifile("unstressed_volumes/vuA", 1.263);

	// modify units of unstressed volume
	Vu[0] = Vu[0]*1e-3;
	Vu[1] = Vu[1]*1e-3;
	Vu[2] = Vu[2]*1e-3;

	// [resistances]
	R[0] = ifile("resistances/rml", 1.021);
	R[1] = ifile("resistances/rlt", 0.3369);
	R[2] = ifile("resistances/rtb", 0.3063);
	R[3] = ifile("resistances/rbA", 0.0817);

	// respiratory params
	resp[0] = ifileC("lungMechanics/rr", 12.);
	resp[1] = ifile("respiratory/ieratio", 0.6);
	resp[2] = ifile("respiratory/frc", 2.4);
	resp[3] = ifileC("pleuralPressure/pl", -5.);
	resp[4] = ifileC("lungMechanics/pmusmin", -5.);
	resp[5] = ifile("respiratory/tauCoeff", 0.2);
	resp[6] = ifile("respiratory/pao", 0.);
	resp[7] = ifile("respiratory/pvent", 0.);
	resp[8] = ifile("respiratory/patm", 0.);
	resp[9] = ifile("respiratory/IAPee", 2.45);	

	if (verbose)
		printLungMechParameters();

}


void lungMechanics::printLungMechParameters() {

	cout << "[timing]" << endl;
	cout << "tIni = " << time << endl;
	// cout << "tEnd = " << timeStop << endl;

	cout << "[compliances]" << endl;
	cout << "cl = " << C[0] << endl;
	cout << "ctr = " << C[1] << endl;
	cout << "cb = " << C[2] << endl;
	cout << "cA = " << C[3] << endl;
	cout << "ccw = " << C[4] << endl;
	cout << "cabd = " << C[5] << endl;


	cout << "[unstressed volumes]" << endl;
	cout << "vul = " << Vu[0] << endl;
	cout << "vutr = " << Vu[1] << endl;
	cout << "vub = " << Vu[2] << endl;
	cout << "vuA = " << Vu[3] << endl;


	cout << "[resistances]" << endl;
	cout << "rml = " << R[0] << endl;
	cout << "rlt = " << R[1] << endl;
	cout << "rtb = " << R[2] << endl;
	cout << "rbA = " << R[3] << endl;


	cout << "[respiratory params]" << endl;
	cout << "rr = " << resp[0] << endl;
	cout << "ieratio = " << resp[1] << endl;
	cout << "frc = " << resp[2] << endl;
	cout << "pplee = " << resp[3] << endl;
	cout << "pmusmin = " << resp[4] << endl;
	cout << "tauCoeff = " << resp[5] << endl;
	cout << "pao = " << resp[6] << endl;
	cout << "pvent = " << resp[7] << endl;
	cout << "patm = " << resp[8] << endl;
	cout << "IAPee= "  << resp[9] << endl;
	
}

void lungMechanics::InitialCond() {
	// [initialConditions]
	vol[0] = Vu[0];
	vol[1] = - C[1] * resp[3] + Vu[1];
	vol[2] = - C[2] * resp[3] + Vu[2];
	vol[3] = - C[3] * resp[3] + Vu[3];
	vol[4] = 0.;
	vol[5] = 0.; 

	if (verbose){
		cout << " -- lungMechanics::saveInitialCond()"<< endl;
		cout << "[initial conditions]" << endl;
		cout << "vl = " << vol[0] << endl;
		cout << "vtr = " << vol[1] << endl;
		cout << "vb = " << vol[2] << endl;
		cout << "vA = " << vol[3] << endl;
		cout << "vpl = " << vol[4] << endl;
		cout << "xi = " << vol[5] << endl;
	}

	getAlgebraicRelations();
}


void lungMechanics::getTimeConstants() {
	T = 60 / resp[0];
	te = T / (resp[1] + 1);
	ti = resp[1] * te;
	tau = resp[5] * te;
}


void lungMechanics::getPmus() {
	getTimeConstants();
	u = vol[5] - floor(vol[5]);

	if (u <= ti / T ){
		Pmus = -1 * resp[4] * pow((u*T), 2) / (te * ti)  + (resp[4] * T * (u*T)) / (ti * te) ; 
	}else{
		Pmus = resp[4] / (1 - exp(- te / tau)) * (exp(-((u*T) - ti) / tau) - exp(-te / tau));
	}
	
	
}


void lungMechanics::getAlgebraicRelations() {
	vd = vol[0] + vol[1] + vol[2];
	v = vd + vol[3];
	// [Pressures] NB: the unit is cmH2O
	getPmus();
	P[4] = 1 / C[4] * (vol[4]) + Pmus + resp[3]; 	// ppl1
	P[0] = (vol[0] - Vu[0]) / C[0]; 				// pl 
	P[1] = (vol[1] - Vu[1]) / C[1] + P[4];			// ptr
	P[2] = (vol[2] - Vu[2]) / C[2] + P[4];			// pb
	P[3] = (vol[3] - Vu[3]) / C[3] + P[4];			// pA
	P[5] = 1 / C[5] * (-vol[4]) - Pmus + resp[9]; 	// abdominal pressure
	Vdot = (resp[6] - P[0]) / R[0];
	VAdot = (P[2] - P[3]) / R[3];
}


void lungMechanics::getTimeDerivative() {
	getAlgebraicRelations();

	dvdt[0] = (resp[6] - P[0]) / R[0] - (P[0] - P[1]) / R[1];	// dvldt
	dvdt[1] = (P[0] - P[1]) / R[1] - (P[1] - P[2]) / R[2];		// dvtrdt
	dvdt[2] = (P[1] - P[2]) / R[2] - (P[2] - P[3]) / R[3];		// dvbdt
	dvdt[3] = (P[2] - P[3]) / R[3];								// dvAdt
	dvdt[4] = (P[0] - P[1]) / R[1] ;							// dvpldt
	dvdt[5] = 1/T; // dxidt
}


void lungMechanics::output() {
	// State Vars
	sampleLungVol << scientific << setprecision(12) << time << " "
			<< vol[0] << " " 
			<< vol[1] << " "
			<< vol[2] << " "
			<< vol[3] << " "
			<< vol[4] << " "
			<< vol[5] << " "
			<< vd << " "
			<< v << " "
			<< "\n";
	sampleLungVol.flush();		

	// Algebraic relations
	sampleLungPr << scientific << setprecision(12) << time << " "
			<< P[0] << " " 
			<< P[1] << " "
			<< P[2] << " "
			<< P[3] << " "
			<< P[4] << " "
			<< Pmus << " "
			<< P[5] << " "
			<< "\n";
	sampleLungPr.flush();
	
	sampleLungFlow << scientific << setprecision(12) << time << " "
			<< Vdot << " " 
			<< VAdot << " "
			<< "\n";
	sampleLungFlow.flush();	
}
