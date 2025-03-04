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

#include "cardiovascular.h"

using namespace std;
using std::scientific;

void cardiovascular::init(string ifile, string ifileC, string outDir) {

	// read parameter file
	if (verbose)
		cout << "Reading CVS params ..." << endl;
	if (ifile != "NOFILE") {
        if (ifileC !="NOFILE"){
		    readCVSparameters(ifile,ifileC);
        }
	}
	if (verbose)
		cout << "Done reading CVS params" << endl;
	
	// intialiaze output file
	nameCVS = "CVS";
    string statenameStateVars = outDir + nameCVS + "_stateVars.out";
	sampleCVSstateVars.open(statenameStateVars.c_str());
    sampleCVSstateVars << "# [0]: time;" << " "
					<< "[1]: vsa;" << " "
                    << "[2]: vep;" << " "
                    << "[3]: vsp;" << " "
                    << "[4]: vmp;" << " "
                    << "[5]: vhp;" << " "
                    << "[6]: vbp;" << " "
                    << "[7]: vev;" << " "
                    << "[8]: vsv;" << " "
                    << "[9]: vmv;" << " "
                    << "[10]: vhv;" << " "
                    << "[11]: vbv;" << " "
                    << "[12]: vtv;" << " "
                    << "[13]: vpa;" << " " 
                    << "[14]: vpp;" << " " 
                    << "[15]: vpv;" << " " 
                    << "[16]: vps;" << " " 
					<< "[17]: vla;" << " " 
					<< "[18]: vlv;" << " " 
                    << "[19]: epsilon;" << " " 
                    << "[20]: vra;" << " " 
                    << "[21]: vrv;" << " " 
				    << "[22]: qTV;" << " " 
                    << "[23]: xiTV;" << " " 
                    << "[24]: qPV;" << " " 
                    << "[25]: xiPV;" << " " 
                    << "[26]: qMV;" << " " 
                    << "[27]: xiMV;" << " " 
                    << "[28]: qAV;" << " " 
                    << "[29]: xiAV;" << " " 
                    << "\n";

	string statenameSys = outDir + nameCVS + "_SysCirc.out";
	sampleCVSsysCirc.open(statenameSys.c_str());
	sampleCVSsysCirc << "# [0]: time;" << " "
                    << "[1]: psa;" << " "
                    << "[2]: qsa;" << " "
                    << "[3]: pep;" << " "
                    << "[4]: psp;" << " "
                    << "[5]: pmp;" << " "
                    << "[6]: php;" << " "
                    << "[7]: pbp;" << " "
                    << "[8]: qep;" << " "
                    << "[9]: qsp;" << " "
                    << "[10]: qmp;" << " "
                    << "[11]: qhp;" << " "
                    << "[12]: qbp;" << " "
                    << "[13]: ptvtm;" << " "
                    << "[14]: ptv;" << " "
                    << "[15]: rtv;" << " "
                    << "[16]: pev;" << " " 
                    << "[17]: psv;" << " "
                    << "[18]: pmv;" << " "
                    << "[19]: phv;" << " "
                    << "[20]: pbv;" << " "
                    << "[21]: qev;" << " "
                    << "[22]: qsv;" << " "
                    << "[23]: qmv;" << " "
                    << "[24]: qhv;" << " "
                    << "[25]: qbv;" << " "
                    << "[26]: qtv;" << " "
                    << "[27]: vsv;" << " "
                    << "[28]: vmv;" << " "                    
                    << "[29]: qepin;" << " "
                    << "[30]: qspin;" << " "
                    << "[31]: qmpin;" << " "
                    << "[32]: qhpin;" << " "
                    << "[33]: qbpin;" << " "
					<< "\n";

	string statenamePul = outDir + nameCVS + "_PulCirc.out";
	sampleCVSpulCirc.open(statenamePul.c_str());
    sampleCVSpulCirc << "# [0]: time;" << " "
					<< "[1]: ppa;" << " " 
					<< "[2]: ppv;" << " "
					<< "[3]: qpp;" << " "
					<< "[4]: qps;" << " "
					<< "[5]: qpv;" << " "
					<< "[6]: qpa;" << " "
                    << "[7]: ppp;" << " "
                    << "[8]: pps;" << " "
                    << "[9]: qppin;" << " "
                    << "[10]: qpsin;" << " "
					<< "\n";

    string statenameValve = outDir + nameCVS + "_Valve.out";
	sampleCVSvalve.open(statenameValve.c_str());
	sampleCVSvalve << "# [0]: time;" << " "
                    << "[1]: ppprox;" << " " 
                    << "[2]: psprox;" << " " 
                    << "[3]: deltapTV;" << " " 
                    << "[4]: deltapPV;" << " " 
                    << "[5]: deltapMV;" << " " 
                    << "[6]: deltapAV;" << " " 
					<< "\n";

    string statenameHeart = outDir + nameCVS + "_Heart.out";
	sampleCVSheart.open(statenameHeart.c_str());
	sampleCVSheart << "# [0]: time;" << " "
                    << "[1]: tsys;" << " " 
                    << "[2]: u;" << " " 
                    << "[3]: phi;" << " " 
                    << "[4]: pmaxlv;" << " " 
                    << "[5]: pmaxrv;" << " " 
                    << "[6]: plv;" << " " 
                    << "[7]: prv;" << " " 
                    << "[8]: pra;" << " " 
                    << "[9]: pla;" << " " 
                    << "[10]: eLA;" << " " 
                    << "[11]: eRA;" << " " 
                    << "[12]: rlv;" << " " 
                    << "[13]: rrv;" << " " 
					<< "\n";

    InitialCond(ifile);
	
	if (verbose){
		cout << "cardiovascular::init :: DONE " << endl;
	}

}


void cardiovascular::readCVSparameters(string _file, string _fileC) {
	if (verbose)
		cout << " -- cardiovascular::readCVSparameters() reading from "
				<< _file << " and "
                << _fileC << endl;

	GetPot ifile(_file.c_str());
	GetPot ifileC(_fileC.c_str());  // common parameters are assignes in commonParams.dat

	// Default values from Albanese et al (2015)
	// if a variable is not found under [category]
	// then a default value is assigned
	
	// [compliances]
	Cparam[0] = ifile("compliances/csa",0.);
	Cparam[1] = ifile("compliances/cep",0.);
	Cparam[2] = ifile("compliances/csp",0.);
	Cparam[3] = ifile("compliances/cmp",0.);
	Cparam[4] = ifile("compliances/chp",0.);
	Cparam[5] = ifile("compliances/cbp",0.);
    Cparam[6] = ifile("compliances/cev",0.);
    Cparam[7] = ifile("compliances/csv",0.);
    Cparam[8] = ifile("compliances/cmv",0.);
    Cparam[9] = ifile("compliances/chv",0.);
    Cparam[10] = ifile("compliances/cbv",0.);
    Cparam[11] = ifile("compliances/cpa",0.);
    Cparam[12] = ifile("compliances/cpp",0.);
    Cparam[13] = ifile("compliances/cpv",0.);

	// [unstressed volumes]
	VuParam[0] = ifile("unstressed_volumes/vusa",0.);
	VuParam[1] = ifile("unstressed_volumes/vuep",0.);
    VuParam[2] = ifile("unstressed_volumes/vusp",0.);
	VuParam[3] = ifile("unstressed_volumes/vump",0.);
    VuParam[4] = ifile("unstressed_volumes/vuhp",0.);
    VuParam[5] = ifile("unstressed_volumes/vubp",0.);
    VuParam[6] = ifileC("cardiovascularControl/vuev0",0.);
    VuParam[7] = ifileC("cardiovascularControl/vusv0",0.);
    VuParam[8] = ifileC("cardiovascularControl/vumv0",0.);
    VuParam[9] = ifile("unstressed_volumes/vuhv",0.);
    VuParam[10] = ifile("unstressed_volumes/vubv",0.);
    VuParam[11] = ifile("unstressed_volumes/vupa",0.);
    VuParam[12] = ifile("unstressed_volumes/vupp",0.);
    VuParam[13] = ifile("unstressed_volumes/vupv",0.);

	// [resistances]
	Rparam[0] = ifile("resistances/Rp",0.);
    Rparam[1] = ifile("resistances/rev",0.);
    Rparam[2] = ifile("resistances/rsv",0.);
    Rparam[3] = ifile("resistances/rmv",0.);
    Rparam[4] = ifile("resistances/rhv",0.);
    Rparam[5] = ifile("resistances/rbv",0.);
    Rparam[6] = ifile("resistances/rpp",0.);
    Rparam[7] = ifile("resistances/rpv",0.);
    Rparam[8] = ifileC("cardiovascularControl/rsp0",0.);
    Rparam[9] = ifileC("cardiovascularControl/rep0",0.);
    Rparam[10] = ifileC("cardiovascularControl/rmp0",0.);
    Rparam[11] = ifileC("cardiovascularControl/rhp0",0.);
    Rparam[12] = ifileC("cardiovascularControl/rbp0",0.);
    Rparam[13] = ifile("resistances/Rpp",0.);
    Rparam[14] = ifile("resistances/Rps",0.);

	// [thoracic veins]
	TveinParam[0] = ifile("thoracicVeins/tvd1",0.);
    TveinParam[1] = ifile("thoracicVeins/tvk1",0.);
    TveinParam[2] = ifile("thoracicVeins/tvkr",0.);
    TveinParam[3] = ifile("thoracicVeins/tvvmax",0.);
    TveinParam[4] = ifile("thoracicVeins/tvvu",0.);
    TveinParam[5] = ifile("thoracicVeins/tvr0",0.);
    TveinParam[6] = ifile("thoracicVeins/tvd2",0.);
    TveinParam[7] = ifile("thoracicVeins/tvk2",0.);
    TveinParam[8] = ifile("thoracicVeins/tvvmin",0.);
    TveinParam[9] = ifile("thoracicVeins/tvkxp",0.);
    TveinParam[10] = ifile("thoracicVeins/tvkxv",0.);


    // [atria]
    atriaParam[0] = ifile("atria/vula",0.); 
    atriaParam[1] = ifile("atria/vura",0.); 
    atriaParam[2] = ifile("atria/s",0.); 
    atriaParam[3] = ifile("atria/trRA",0.); 
    atriaParam[4] = ifile("atria/trLA",0.); 
    atriaParam[5] = ifile("atria/tcRA",0.); 
    atriaParam[6] = ifile("atria/tcLA",0.); 
    atriaParam[7] = ifile("atria/TrpRA",0.); 
    atriaParam[8] = ifile("atria/TrpLA",0.); 
    atriaParam[9] = ifile("atria/TcpRA",0.); 
    atriaParam[10] = ifile("atria/TcpLA",0.); 
    atriaParam[11] = ifile("atria/eALA",0.); 
    atriaParam[12] = ifile("atria/eBLA",0.); 
    atriaParam[13] = ifile("atria/eARA",0.); 
    atriaParam[14] = ifile("atria/eBRA",0.); 

    // [ventricles]
    ventrParam[0] = ifile("ventricles/p0lv",0.);
    ventrParam[1] = ifile("ventricles/p0rv",0.);
    ventrParam[2] = ifile("ventricles/kelv",0.);
    ventrParam[3] = ifile("ventricles/kerv",0.);
    ventrParam[4] = ifile("ventricles/vulv",0.);
    ventrParam[5] = ifile("ventricles/vurv",0.);
    ventrParam[6] = ifileC("cardiovascularControl/emaxlv0",0.);
    ventrParam[7] = ifileC("cardiovascularControl/emaxrv0",0.);
    ventrParam[8] = ifile("ventricles/krlv",0.);
    ventrParam[9] = ifile("ventricles/krrv", 0.0014);

    // [heart valves]
    valveParam[0] = ifile("valve/AmaxTV",-1.);
    valveParam[1] = ifile("valve/AminTV",-1.);
    valveParam[2] = ifile("valve/lTV",0.);
    valveParam[3] = ifile("valve/deltap_openTV",0.);
    valveParam[4] = ifile("valve/k_openTV",0.);
    valveParam[5] = ifile("valve/deltap_closeTV",0.);
    valveParam[6] = ifile("valve/k_closeTV",0.);
    valveParam[7] = ifile("valve/AmaxPV",-1.);
    valveParam[8] = ifile("valve/AminPV",-1.);
    valveParam[9] = ifile("valve/lPV",0.);
    valveParam[10] = ifile("valve/deltap_openPV",0.);
    valveParam[11] = ifile("valve/k_openPV",0.);
    valveParam[12] = ifile("valve/deltap_closePV",0.);
    valveParam[13] = ifile("valve/k_closePV",0.);
    valveParam[14] = ifile("valve/AmaxMV",-1.);
    valveParam[15] = ifile("valve/AminMV",-1.);
    valveParam[16] = ifile("valve/lMV",0.);
    valveParam[17] = ifile("valve/deltap_openMV",0.);
    valveParam[18] = ifile("valve/k_openMV",0.);
    valveParam[19] = ifile("valve/deltap_closeMV",0.);
    valveParam[20] = ifile("valve/k_closeMV",0.);
    valveParam[21] = ifile("valve/AmaxAV",-1.);
    valveParam[22] = ifile("valve/AminAV",-1.);
    valveParam[23] = ifile("valve/lAV",0.);
    valveParam[24] = ifile("valve/deltap_openAV",0.);
    valveParam[25] = ifile("valve/k_openAV",0.);
    valveParam[26] = ifile("valve/deltap_closeAV",0.);
    valveParam[27] = ifile("valve/k_closeAV",0.);
    valveParam[28] = ifile("valve/rpprox",0.);
    valveParam[29] = ifile("valve/rsprox",0.);
    valveParam[30] = ifile("valve/ATV",0.);
    valveParam[31] = ifile("valve/APV",0.);
    valveParam[32] = ifile("valve/AMV",0.);
    valveParam[33] = ifile("valve/AAV",0.);

    // Assign max and min valve areas
    valveParam[0] = valveParam[30];        // AmaxTV
    valveParam[1] = valveParam[30]*1e-5;   // AminTV
    valveParam[7] = valveParam[31];        // AmaxPV
    valveParam[8] = valveParam[31]*1e-5;   // AminPV
    valveParam[14] = valveParam[32];       // AmaxMV
    valveParam[15] = valveParam[32]*1e-5;  // AminMV
    valveParam[21] = valveParam[33];       // AmaxAV
    valveParam[22] = valveParam[33]*1e-5;  // AminAV

    // [additional parameters]
    params[0] = ifileC("shunt/sh",0.); 
    params[1] = ifile("param/rho", 1.04); 
    params[2] = ifile("param/T", 0.833);  
    params[3] = ifile("param/tsys0",0.); 
    params[4] = ifile("param/ksys",0.); 
	params[5] = ifileC("pleuralPressure/pl", 1000.); 

    // controlled variables
    controlledVars[0] = params[2];      // T 
    controlledVars[1] = VuParam[6];     // vuev
    controlledVars[2] = VuParam[7];     // vusv
    controlledVars[3] = VuParam[8];     // vumv
    controlledVars[4] =  Rparam[9];     // rep
    controlledVars[5] = Rparam[8];      // rsp
    controlledVars[6] = Rparam[10];     // rmp
    controlledVars[7] = Rparam[11];     // rhp
    controlledVars[8] = Rparam[12];     // rbp
    controlledVars[9] = ventrParam[6];  // emaxlv
    controlledVars[10] = ventrParam[7]; // emaxrv

	if (verbose)
		printCVSparameters();
}

void cardiovascular::InitialCond(string _file) {
	if (verbose)
		cout << " -- cardiovascular::saveInitialCond() reading from "
				<< _file << endl;

	GetPot ifile(_file.c_str());
	// [initialConditions]
    stateVars[0] = ifile("initial_conditions/vsa",0.);
    stateVars[1] = ifile("initial_conditions/vep",0.);
    stateVars[2] = ifile("initial_conditions/vsp",0.);
    stateVars[3] = ifile("initial_conditions/vmp",0.);
    stateVars[4] = ifile("initial_conditions/vhp",0.);
    stateVars[5] = ifile("initial_conditions/vbp",0.);
    stateVars[6] = ifile("initial_conditions/vev",0.);
    stateVars[7] = ifile("initial_conditions/vsv",0.);
    stateVars[8] = ifile("initial_conditions/vmv",0.);
    stateVars[9] = ifile("initial_conditions/vhv",0.);
    stateVars[10] = ifile("initial_conditions/vbv",0.);
    stateVars[11] = ifile("initial_conditions/vtv",0.);
    stateVars[12] = ifile("initial_conditions/vpa",0.);
    stateVars[13] = ifile("initial_conditions/vpp",0.);
    stateVars[14] = ifile("initial_conditions/vpv",0.);
    stateVars[15] = ifile("initial_conditions/vps",0.);
    stateVars[16] = ifile("initial_conditions/vla",0.);
    stateVars[17] = ifile("initial_conditions/vlv",0.);
    stateVars[18] = ifile("initial_conditions/epsilon",0.);
    stateVars[19] = ifile("initial_conditions/vra",0.);
    stateVars[20] = ifile("initial_conditions/vrv",0.);
    stateVars[21] = ifile("initial_conditions/qTV",0.);
    stateVars[22] = ifile("initial_conditions/xiTV",0.);
    stateVars[23] = ifile("initial_conditions/qPV",0.);
    stateVars[24] = ifile("initial_conditions/xiPV",0.);
    stateVars[25] = ifile("initial_conditions/qMV",0.);
    stateVars[26] = ifile("initial_conditions/xiMV",0.);
    stateVars[27] = ifile("initial_conditions/qAV",0.);
    stateVars[28] = ifile("initial_conditions/xiAV",0.);

    getAlgebraicRelations();
}


void cardiovascular::printCVSparameters() {
    cout << "[controlled variables]" << endl;
    for (int i = 0; i < controlledVars.size(); i++) {
		cout << i << " value: " << controlledVars[i] << endl;        
	}
    cout << "[additional params]" << endl;
    for (int i = 0; i < params.size(); i++) {
        cout << i << " value: " << params[i] << endl;        
	}
    cout << "[valve params]" << endl;
    for (int i = 0; i < valveParam.size(); i++) {
        cout << i << " value: " << valveParam[i] << endl;        
	}
    cout << "[ventricles params]" << endl;
    for (int i = 0; i < ventrParam.size(); i++) {
        cout << i << " value: " << ventrParam[i] << endl;        
	}
    cout << "[atria params]" << endl;
    for (int i = 0; i < atriaParam.size(); i++) {
        cout << i << " value: " << atriaParam[i] << endl;        
	}
    cout << "[thoracic veins params]" << endl;
    for (int i = 0; i < TveinParam.size(); i++) {
        cout << i << " value: " << TveinParam[i] << endl;        
	}
    cout << "[resistances]" << endl;
    for (int i = 0; i < Rparam.size(); i++) {
        cout << i << " value: " << Rparam[i] << endl;        
	}
    cout << "[unstressed volumes]" << endl;
    for (int i = 0; i < VuParam.size(); i++) {
        cout << i << " value: " << VuParam[i] << endl;        
	}
    cout << "[compliances]" << endl;
    for (int i = 0; i < Cparam.size(); i++) {
        cout << i << " value: " << Cparam[i] << endl;        
	}
}

double cardiovascular::getqValvedt(double deltap, double Amax, double Amin, double q, double l, double xi) {
    double A;
    double b; 
    double L;
    double dqdt;
    A = (Amax-Amin)*xi + Amin;         
    b = 0.5*params[1]/pow(A,2.);                // M - Bernoulli resistance 
    L = params[1]*l/A;                          // M - inertance
    dqdt=(deltap*1333.22 -b*q*abs(q))/L; 
    return dqdt;
}

double cardiovascular::getxiValvesdt(double deltap,double deltap_open, double k_open, double deltap_close,double k_close,double xi) {
    double dxidt;                  // M - variation of the opening state of valve
    dxidt=0.;
    if (deltap>=deltap_open){
        dxidt=(1-xi)*k_open*(deltap-deltap_open); 
    }        
    else if (deltap<deltap_close){
        dxidt=xi*k_close*(deltap-deltap_close);
    }
    return dxidt;
}

double cardiovascular::geteAtria(double t, double tr, double Trp,double tc,double Tcp) {
    double e; 
    if (t<=tr+Trp-1.){
        e=0.5*(1.+cos(M_PI*(t+1.-tr)/Trp));
    }
    else if (t<=tc){ 
        e=0.;
    }
    else if (t<=tc+Tcp){ 
        e=0.5*(1.-cos(M_PI*(t-tc)/Tcp));
    }
    else if (t<=1.){
        e=0.5*(1.+cos(M_PI*(t-tr)/Trp));
    }
    return e;
}

void cardiovascular::getAlgebraicRelations() {
    // SYSTEMIC CIRCULATION
    // systemic arteries cyrculation
    sysCircAlg[0] = (stateVars[0]-VuParam[0])/Cparam[0];    // A3 - pressure in systemic arteries
    // pheripheral compartments cyrculation   
    sysCircAlg[2] = (stateVars[1]-VuParam[1])/Cparam[1];    // A5 - pressure e
    sysCircAlg[3] = (stateVars[2]-VuParam[2])/Cparam[2];    // A5 - pressure s
    sysCircAlg[4] = (stateVars[3]-VuParam[3])/Cparam[3];    // A5 - pressure m
    sysCircAlg[5] = (stateVars[4]-VuParam[4])/Cparam[4];    // A5 - pressure h
    sysCircAlg[6] = (stateVars[5]-VuParam[5])/Cparam[5];    // A5 - pressure b

    sysCircAlg[28] = ((sysCircAlg[0]-sysCircAlg[2])/Rparam[0]);     // our modification - flow e
    sysCircAlg[29] = ((sysCircAlg[0]-sysCircAlg[3])/Rparam[0]);     // our modification - flow s
    sysCircAlg[30] = ((sysCircAlg[0]-sysCircAlg[4])/Rparam[0]);     // our modification - flow m
    sysCircAlg[31] = ((sysCircAlg[0]-sysCircAlg[5])/Rparam[0]);     // our modification - flow h
    sysCircAlg[32] = ((sysCircAlg[0]-sysCircAlg[6])/Rparam[0]);     // our modification - flow b
    // systemic arteries cyrculation
    sysCircAlg[1] = sysCircAlg[28] + sysCircAlg[29] + sysCircAlg[30] + sysCircAlg[31] + sysCircAlg[32];     // our modification - flow in systemic arteries
    // systemic venous circulation
    sysCircAlg[26] = stateVars[7];     // A14 - volume s
    sysCircAlg[27] = stateVars[8];     // A14 - volume m
    sysCircAlg[15] = max(stateVars[6]-controlledVars[1],0.)/Cparam[6];  // A - pressure s
    sysCircAlg[16] = (sysCircAlg[26]-controlledVars[2])/Cparam[7];      // A - pressure e
    sysCircAlg[17] = (sysCircAlg[27]-controlledVars[3])/Cparam[8];      // A - pressure m
    sysCircAlg[18] = (stateVars[9]-VuParam[9])/Cparam[9];               // A - pressure h
    sysCircAlg[19] = (stateVars[10]-VuParam[10])/Cparam[10];            // A - pressure b
    // pheripheral compartments cyrculation
    sysCircAlg[7] = max(sysCircAlg[2]-sysCircAlg[15],0.)/controlledVars[4];     // A6 - flow e
    sysCircAlg[8] = max(sysCircAlg[3]-sysCircAlg[16],0.)/controlledVars[5];     // A6 - flow s
    sysCircAlg[9] = max(sysCircAlg[4]-sysCircAlg[17],0.)/controlledVars[6];     // A6 - flow m 
    sysCircAlg[10] = max(sysCircAlg[5]-sysCircAlg[18],0.)/controlledVars[7];    // A6 - flow h
    sysCircAlg[11] = max(sysCircAlg[6]-sysCircAlg[19],0.)/controlledVars[8];    // A6 - flow b
    // thoracic veins
    double psi;
    psi = TveinParam[9]/exp(stateVars[11]/TveinParam[10]-1.) ;
    if ((stateVars[11])>=TveinParam[4]){
        sysCircAlg[12] = TveinParam[0] + TveinParam[1]*(stateVars[11]-TveinParam[4]) - psi;
    }
    else{
        sysCircAlg[12] = TveinParam[6] + TveinParam[7]*exp(stateVars[11]/TveinParam[8])-psi;
    }
    sysCircAlg[13] = sysCircAlg[12] + params[5];                                                                    // A17 - pressure 
    sysCircAlg[14] = TveinParam[2]*pow(TveinParam[3]/stateVars[11],2)+TveinParam[5];                // resistance
    sysCircAlg[25] = max(0.,(sysCircAlg[13]-heartAlg[7])/sysCircAlg[14]);                                           // A16 - flow
    // systemic venous circulation
    sysCircAlg[20] = max(sysCircAlg[15]-sysCircAlg[13],0.)/Rparam[1];       // A13 - flow e
    sysCircAlg[21] = max(sysCircAlg[16]-sysCircAlg[13],0.)/Rparam[2];       // A13 - flow s
    sysCircAlg[22] = max((sysCircAlg[17]-sysCircAlg[13])/Rparam[3],0.);     // A13 - flow m
    sysCircAlg[23] = max(sysCircAlg[18]-sysCircAlg[13],0.)/Rparam[4];       // A13 - flow h
    sysCircAlg[24] = max((sysCircAlg[19]-sysCircAlg[13])/Rparam[5],0.);     // A13 - flow b

    // PULMONARY CIRCULATION    
    pulCircAlg[0] = (stateVars[12] -VuParam[11] )/Cparam[11]+ params[5];                                            // A20 - pressure arteries
    pulCircAlg[6] = (stateVars[13] - VuParam[12]*((100-params[0])/100))*100/Cparam[12]/(100-params[0])+params[5];   // A - pressure peripheral
    pulCircAlg[7] = (stateVars[15] - VuParam[12]*((params[0])/100))*100/Cparam[12]/params[0]+params[5];             // A - pressure shunt
    pulCircAlg[1] = (stateVars[14]-VuParam[13])/Cparam[13] + params[5];                                             // A29 - pressure veins
    pulCircAlg[8] = (pulCircAlg[0]-pulCircAlg[6])/Rparam[13];                                                       // our modification - flow in peripheral
    pulCircAlg[9] = (pulCircAlg[0]-pulCircAlg[7])/Rparam[14];                                                       // our modification - flow in shunt
    pulCircAlg[5] = pulCircAlg[8] + pulCircAlg[9];                                                                  // flow through arteries
    pulCircAlg[2] = (pulCircAlg[6]-pulCircAlg[1])/(Rparam[6]*100/(100-params[0]));                                  // A22 - flow through peripheral
    pulCircAlg[3] = (pulCircAlg[7]-pulCircAlg[1])/(Rparam[6]*100/params[0]);                                        // A23 - flow through shunt
    pulCircAlg[4] = max((pulCircAlg[1]-heartAlg[8]),0.)/Rparam[7];   
    
    // HEART
	// general heart
    heartAlg[0] = params[3]-params[4]/controlledVars[0];    // U22 - systolic period
    heartAlg[1] = stateVars[18]-floor(stateVars[18]);       // U21 - fraction of cardiac cycle u
    heartAlg[2] = 0.;                                       // U19 - phi
    if (heartAlg[1] <= heartAlg[0]/controlledVars[0]){
    heartAlg[2] = pow(sin(M_PI*controlledVars[0]/heartAlg[0]*heartAlg[1]),2);
    }
    // ventricles
    heartAlg[3] = heartAlg[2]*controlledVars[9]*(stateVars[17]-ventrParam[4]) + (1.-heartAlg[2])*ventrParam[0]*(exp(ventrParam[2]*stateVars[17])-1.) + params[5]; // U18 - isovolumetric pressure LV
    heartAlg[11] = max(0.0001,ventrParam[8]*heartAlg[3]);   // U16 - resistance LV
    heartAlg[5] = heartAlg[3]-heartAlg[11]*stateVars[27];   // U17 - pressure LV
    
    heartAlg[4] = heartAlg[2]*controlledVars[10]*(stateVars[20]-ventrParam[5]) + (1.-heartAlg[2])*ventrParam[1]*(exp(ventrParam[3]*stateVars[20])-1.) + params[5]; // U29 - isovolumetric pressure RV
	heartAlg[12] = max(0.0001,ventrParam[9]*heartAlg[4]);   // U27 - resistance RV
    heartAlg[6] = heartAlg[4]-heartAlg[12]*stateVars[23];   // U28 - pressure RV
    // atria 
    heartAlg[9]=geteAtria(heartAlg[1],atriaParam[4],atriaParam[8],atriaParam[6],atriaParam[10]);    // L - elastance LA
    heartAlg[10]=geteAtria(heartAlg[1],atriaParam[3],atriaParam[7],atriaParam[5],atriaParam[9]);    // L - elastance RA

    heartAlg[7]=(atriaParam[13]*heartAlg[10]+atriaParam[14])*(stateVars[19]-atriaParam[1])+atriaParam[2]*(sysCircAlg[25]-stateVars[21])+params[5];  // L - pressure RA
    heartAlg[8]=(atriaParam[11]*heartAlg[9]+atriaParam[12])*(stateVars[16]-atriaParam[0])+atriaParam[2]*(pulCircAlg[4]-stateVars[25])+params[5];    // L - pressure LA
    // heart valves
    valveAlg[0]=pulCircAlg[0]+valveParam[28]*stateVars[23];   // M - proximal pulmonary pressure
    valveAlg[1]=sysCircAlg[0]+valveParam[29]*stateVars[27];   // M - proximal systemic pressure
    valveAlg[2]=heartAlg[7]-heartAlg[6];                      // M - delta pressure TV 
    valveAlg[3]=heartAlg[6]-valveAlg[0];                      // M - delta pressure PV
    valveAlg[4]=heartAlg[8]-heartAlg[5];                      // M - delta pressure MV
    valveAlg[5]=heartAlg[5]-valveAlg[1];                      // M - delta pressure AV 
}

void cardiovascular::getTimeDerivative() {
	// definition of the ODE system to be solved
    getAlgebraicRelations();
    // systemic arteries circulation
    dvdt[0] = stateVars[27] - sysCircAlg[1];          // A1 - volume

    // systemic peripheral circulation
    dvdt[1] = sysCircAlg[28] - sysCircAlg[7];        // A4 - volume e
    dvdt[2] = sysCircAlg[29] - sysCircAlg[8];        // A4 - volume s
    dvdt[3] = sysCircAlg[30] - sysCircAlg[9];        // A4 - volume m
    dvdt[4] = sysCircAlg[31] - sysCircAlg[10];       // A4 - volume h
    dvdt[5] = sysCircAlg[32] - sysCircAlg[11];       // A4 - volume b

    // systemic venous circulation
    dvdt[6] = sysCircAlg[7] - sysCircAlg[20];         // A - volume e 
    dvdt[7] = sysCircAlg[8] - sysCircAlg[21];         // A - volume s
    dvdt[8] = sysCircAlg[9] - sysCircAlg[22];         // A - volume m
    dvdt[9] = sysCircAlg[10] - sysCircAlg[23];        // A10 - volume h 
    dvdt[10] = sysCircAlg[11] - sysCircAlg[24];       // A11 - volume b
    
    // thoracic veins
    dvdt[11] = sysCircAlg[20] + sysCircAlg[21] + sysCircAlg[22] + sysCircAlg[23] + sysCircAlg[24] - sysCircAlg[25]; // A15 - volume
        
    // pulmonary circulation
    dvdt[12] = stateVars[23] - pulCircAlg[5];                 // A18 - volume arteries
    dvdt[13] = pulCircAlg[8] - pulCircAlg[2];                 // A21 - volume peripheral
    dvdt[14] = pulCircAlg[2] + pulCircAlg[3] - pulCircAlg[4]; // A27 - volume veins
    dvdt[15] = pulCircAlg[9] - pulCircAlg[3];                 // A - volume shunt

    // heart
    dvdt[16] = pulCircAlg[4]-stateVars[25];    // U - volume LA
    dvdt[17] = stateVars[25]-stateVars[27];    // U - volume LV
    dvdt[18]  = 1./controlledVars[0];          // U21 - epsilon
    dvdt[19] = sysCircAlg[25]-stateVars[21];   // U - volume RA
    dvdt[20] = stateVars[21]-stateVars[23];    // U - volume RV

    // valves
    dvdt[21] =getqValvedt(valveAlg[2],valveParam[0],valveParam[1],stateVars[21],valveParam[2],stateVars[22]);        // M - flow through TV
    dvdt[22] =getxiValvesdt(valveAlg[2],valveParam[3],valveParam[4],valveParam[5],valveParam[6],stateVars[22]);      // M - openning state xi TV

    dvdt[23] =getqValvedt(valveAlg[3],valveParam[7],valveParam[8],stateVars[23],valveParam[9],stateVars[24]);        // M - flow through PV
    dvdt[24] =getxiValvesdt(valveAlg[3],valveParam[10],valveParam[11],valveParam[12],valveParam[13],stateVars[24]);  // M - openning state xi PV

    dvdt[25] =getqValvedt(valveAlg[4],valveParam[14],valveParam[15],stateVars[25],valveParam[16],stateVars[26]);     // M - flow through MV
    dvdt[26] =getxiValvesdt(valveAlg[4],valveParam[17],valveParam[18],valveParam[19],valveParam[20],stateVars[26]);  // M - openning state xi MV
    
    dvdt[27] =getqValvedt(valveAlg[5],valveParam[21],valveParam[22],stateVars[27],valveParam[23],stateVars[28]);     // M - flow through AV
    dvdt[28] =getxiValvesdt(valveAlg[5],valveParam[24],valveParam[25],valveParam[26],valveParam[27],stateVars[28]);  // M - openning state xi AV
}

void cardiovascular::output(){
    // write output
    sampleCVSstateVars << scientific << setprecision(12)   ;
    sampleCVSstateVars << time << " ";
    // "state" variables
    for (int i = 0; i < stateVars.size(); i++) {
        sampleCVSstateVars << stateVars[i] << " ";
    }
    sampleCVSstateVars << "\n";
    sampleCVSstateVars.flush();
    // algebraic variables
    sampleCVSsysCirc << scientific << setprecision(12)   ;
    sampleCVSsysCirc << time << " ";
    
    for (int i = 0; i < sysCircAlg.size(); i++) {
        sampleCVSsysCirc << sysCircAlg[i] << " ";
    }
    sampleCVSsysCirc << "\n";
    sampleCVSsysCirc.flush();
    sampleCVSpulCirc << scientific << setprecision(12)   ;
    sampleCVSpulCirc << time << " ";
    // algebraic variables
    for (int i = 0; i < pulCircAlg.size(); i++) {
        sampleCVSpulCirc << pulCircAlg[i] << " ";
    }
    sampleCVSpulCirc << "\n";
    sampleCVSpulCirc.flush();
    sampleCVSvalve << scientific << setprecision(12)   ;
    sampleCVSvalve << time << " ";
    // algebraic variables
    for (int i = 0; i < valveAlg.size(); i++) {
        sampleCVSvalve << valveAlg[i] << " ";
    }
    sampleCVSvalve << "\n";
    sampleCVSvalve.flush();
    sampleCVSheart << scientific << setprecision(12)   ;
    sampleCVSheart << time << " ";
    // algebraic variables
    for (int i = 0; i < heartAlg.size(); i++) {
        sampleCVSheart << heartAlg[i] << " ";
    }
    sampleCVSheart << "\n";
    sampleCVSheart.flush();
}