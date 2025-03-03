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

#include "respiratoryControl.h"

using namespace std;
using std::scientific;

void respiratoryControl::init(string ifile, string ifileC, string outDir){
    if (verbose)
		cout << "Reading Respiratory control params ..." << endl;
	if (ifile != "NOFILE") {
        if (ifileC !="NOFILE"){
		    readRespControlParameters(ifile,ifileC);
        }
	}
	if (verbose)
		cout << "Done reading respiratory control params" << endl;

    nameControl = "RespControl";

    string statenameVol1 = outDir + nameControl + "_algebraicVars.txt";
    sampleAlgVars.open(statenameVol1.c_str());
    sampleAlgVars << "# [0] : time:" << " "
                    << "[1] : Pmusmin " << " "
                    << "[2] : RR " << " "
                    << "[3] : uc " << " "
                    << "[4] : up " << " " 
                    << "\n";


    string statenameVol2 = outDir + nameControl + "_stateVars.txt";
    sampleStateVars.open(statenameVol2.c_str());
    sampleStateVars << "# [0] : time:" << " "
                    << "[1] : DeltaPmusminc " << " "
                    << "[2] : DeltaPmusminp " << " "
                    << "[3] : DeltaRRc " << " "
                    << "[4] : DeltaRRp " << " "
                    << "\n";
	
    InitialCond(ifile);

	if (verbose){
		cout << "RespControlAlbanese::init :: DONE " << endl;
	}
}


void respiratoryControl::readRespControlParameters(string _file, string _fileC){
    if (verbose)
		cout << " -- respiratoryControl::readCardioControlparameters() reading from "
				<< _file << " and "
                << _fileC << endl;

	GetPot ifile(_file.c_str());
	GetPot ifileC(_fileC.c_str());  // common parameters are assignes in commonParams.dat

    Params[0] = ifile("peripheral/Dp", 7.);
    Params[1] = ifile("peripheral/Gpa", 1200.);
    Params[2] = ifile("peripheral/Gpf", 0.8913);
    Params[3] = ifile("peripheral/taupa1", 83.);
    Params[4] = ifile("peripheral/taupa2", 10.);
    Params[5] = ifile("peripheral/taupf1", 174.78);
    Params[6] = ifile("peripheral/taupf2", 17.526);
    Params[7] = ifile("peripheral/fapcn", 3.7);
    Params[8] = ifile("central/Dc", 8.);
    Params[9] = ifile("central/Gca1", 850.);
    Params[10] = ifile("central/Gca2", 0.);
    Params[11] = ifile("central/Gcf1", 0.9);
    Params[12] = ifile("central(Gcf2", 0.);
    Params[13] = ifile("central/tauca1", 105.);
    Params[14] = ifile("central/tauca2", 30.);
    Params[15] = ifile("central/taucf1", 400.);
    Params[16] = ifile("central/taucf2", 35.);
    Params[17] = ifileC("cardiorespiratoryControl/paco2n", 40.);
    Params[18] = ifileC("lungMechanics/pmusmin", -5.);
    Params[19] = ifileC("lungMechanics/rr", 12.);

    if (verbose)
		printRespControlParameters();
}

void respiratoryControl::InitialCond(string _file){

    GetPot ifile(_file.c_str());

    // [initial conditions]
    stateVars[0] = ifile("ics/DeltaPmusminp", 0.);
    stateVars[1] = ifile("ics/DeltaPmusminc", 0.);
    stateVars[2] = ifile("ics/DeltaRRp", 0.);
    stateVars[3] = ifile("ics/DeltaRRc", 0.); 

    getAlgebraicRelations();
}

void respiratoryControl::printRespControlParameters(){
    cout << "[peripheral]"<<endl;
    for (int i = 0; i < 8; i++){
        cout << i << " value: " << Params[i] << endl;
    }

    cout << "[central]"<<endl;
    for (int i = 8; i < 17; i++){
        cout << i << " value: " << Params[i] << endl;
    }

    cout << "[other]"<<endl;
    for (int i = 17; i < 20; i++){
        cout << i << " value: " << Params[i] << endl;
    }
}


void respiratoryControl::getAlgebraicRelations(){
    algVars[0] = - max(0., - Params[18] + (stateVars[0] + stateVars[1])/(980.));
    algVars[1] = Params[19] + stateVars[2] + stateVars[3];
    algVars[2] = (delayedVars[0] - Params[17]);
    algVars[3] = (delayedVars[1] - Params[7]);
}

void respiratoryControl::getTimeDerivative(){
    getAlgebraicRelations();

    double taupa;
    double taupf;
    if (algVars[3] >= 0.){
        taupa = Params[3];
        taupf = Params[5];
    }else{
        taupa = Params[4];
        taupf = Params[6];
    }

    double Gca;
    double Gcf;
    double tauca;
    double taucf;
    if (algVars[2] > 5.){
        Gca = Params[9];
        Gcf = Params[11];
        tauca = Params[13];
        taucf = Params[15];
    }else{
        Gca = Params[10];
        Gcf = Params[12];
        tauca = Params[14];
        taucf = Params[16];

    }

    dvdt[0] = (- stateVars[0] + Params[1] * max(algVars[3],0.)) / taupa;
    dvdt[1] = (- stateVars[1] + Gca * algVars[2]) / tauca;
    dvdt[2] = (- stateVars[2] + Params[2] * max(algVars[3],0.)) / taupf;
    dvdt[3] = (- stateVars[3] + Gcf * algVars[2]) / taucf;
}

void respiratoryControl::output(){
    sampleStateVars << scientific << setprecision(12) << time << " "
                    << stateVars[0] << " "
                    << stateVars[1] << " "
                    << stateVars[2] << " "
                    << stateVars[3] << " "
                    << "\n";    
    sampleStateVars.flush();

    sampleAlgVars << scientific << setprecision(12) << time << " "
                    << algVars[0] << " "
                    << algVars[1] << " "
                    << algVars[2] << " "
                    << algVars[3] << " "
                    << "\n";
    sampleAlgVars.flush();
}

