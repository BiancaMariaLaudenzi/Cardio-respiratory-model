
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

#include "cardiovascularControl.h"

using namespace std;
using std::scientific;
// using std::max;

void cardiovascularControl::init(string ifile, string ifileC, string outDir){
    if (verbose)
		cout << "Reading cardiovascular control params ..." << endl;
	if (ifile != "NOFILE") {
        if (ifileC !="NOFILE"){
		    readCardioControlParameters(ifile,ifileC);
        }
	}
	if (verbose)
		cout << "Done reading cardiovascular control params" << endl;

    // intialiaze output file
    nameControl = "CardioControl";

    string statenameVol1 = outDir + nameControl + "_algebraicVars.txt";
    sampleAlgVars.open(statenameVol1.c_str());
    sampleAlgVars << "# [0]: time;" << " " 
                    << "[1]: gbp;" << " "
                    << "[2]: phib;" << " "
                    << "[3]: phih;" << " "
                    << "[4]: phim;" << " "
                    << "[5]: rbp;" << " "
                    << "[6]: rhp;" << " "
                    << "[7]: rmp;" << " "
                    << "[8]: fab;" << " "
                    << "[9]: Xo2;" << " "
                    << "[10]: phistat;" << " "
                    << "[11]: fcstat;" << " "
                    << "[12]: fcdyn;" << " "
                    << "[13]: fapc;" << " "
                    << "[14]: omegash;" << " "
                    << "[15]: omegasp;" << " "
                    << "[16]: omegasv;" << " "
                    << "[17]: thetash;" << " "
                    << "[18]: thetasp;" << " "
                    << "[19]: thetasv;" << " "
                    << "[20]: phiap;" << " "
                    << "[21]: fsh;" << " "
                    << "[22]: fsp;" << " "
                    << "[23]: fsv;" << " "
                    << "[24]: fv;" << " "
                    << "[25]: sigmarmp;" << " "
                    << "[26]: sigmarsp;" << " "
                    << "[27]: sigma_rep;" << " "
                    << "[28]: sigmavumv;" << " "
                    << "[29]: sigmavusv;" << " "
                    << "[30]: sigmavuev;" << " "
                    << "[31]: sigmaemaxlv;" << " "
                    << "[32]: sigmaemaxrv;" << " "
                    << "[33]: sigmaTs;" << " "
                    << "[34]: sigmaTv;" << " "
                    << "[35]: rmp;" << " "
                    << "[36]: rsp;" << " "
                    << "[37]: rep;" << " "
                    << "[38]: vumv;" << " "
                    << "[39]: vusv;" << " "
                    << "[40]: vuev;" << " "
                    << "[41]: emaxlv;" << " "
                    << "[42]: emaxrv;" << " "
                    << "[43]: T;" << " "
                    << "[44]: rbp;" << " "
                    << "[45]: rhp;" << " "
                    << "\n";


    string statenameVol2 = outDir + nameControl + "_differentialVars.txt";
    sampleDiffVars.open(statenameVol2.c_str());
    sampleDiffVars << "# [0]: time;" << " "
                      << "[1]: xbO2;" << " "
                      << "[2]: xbCO2;" << " "
                      << "[3]: xhO2;" << " "
                      << "[4]: xhCO2;" << " "
                      << "[5]: xmO2;" << " "
                      << "[6]: xmCO2;" << " "
                      << "[7]: Ptilde;" << " "
                      << "[8]: phico2Dyn;" << " "
                      << "[9]: phiapc;" << " "
                      << "[10]: fap;" << " "
                      << "[11]: deltathetao2sh;" << " "
                      << "[12]: deltathetaco2sh;" << " "
                      << "[13]: deltathetao2sp;" << " "
                      << "[14]: deltathetaco2sp;" << " "
                      << "[15]: deltathetao2sv;" << " "
                      << "[16]: deltathetaco2sv;" << " "
                      << "[17]: deltarmp;" << " "
                      << "[18]: deltarsp;" << " "
                      << "[19]: deltarep;" << " "
                      << "[20]: deltavumv;" << " "
                      << "[21]: deltavusv;" << " "
                      << "[22]: deltavuev;" << " "
                      << "[23]: deltaemaxlv;" << " "
                      << "[24]: deltaemaxrv;" << " "
                      << "[25]: deltaTs;" << " "
                      << "[26]: deltaTv;" << " "
                      << "\n";
    
    InitialCond(ifile);

	if (verbose){
		cout << "CardioControlAlbanese::init :: DONE " << endl;
	}
}


void cardiovascularControl::readCardioControlParameters(string _file, string _fileC){
    if (verbose)
		cout << " -- cardiovascularControl::readCardioControlparameters() reading from "
				<< _file << " and "
                << _fileC << endl;

    
	GetPot ifile(_file.c_str());
	GetPot ifileC(_fileC.c_str());  // common parameters are assignes in commonParams.dat

    // basal Heart Rate
    T0 = ifile("basalHR/T0", 0.58);
    
    // [commonParams]
    commonParams[0] = ifileC("cardiovascularControl/rmp0", 2.106);
    commonParams[1] = ifileC("cardiovascularControl/rsp0", 1.8675);
    commonParams[2] = ifileC("cardiovascularControl/rep0", 1.24125);    
    commonParams[3] = ifileC("cardiovascularControl/vumv0", 503.26);
    commonParams[4] = ifileC("cardiovascularControl/vusv0", 1435.4);
    commonParams[5] = ifileC("cardiovascularControl/vuev0", 640.73);
    commonParams[6] = ifileC("cardiovascularControl/emaxlv0", 2.3920);
    commonParams[7] = ifileC("cardiovascularControl/emaxrv0", 1.412);
    commonParams[8] = ifileC("cardiovascularControl/rbp0", 6.6667);
    commonParams[9] = ifileC("cardiovascularControl/rhp0", 19.71);

    // [commonControlParams]
    commonControlParams[0] = ifileC("cardiorespiratoryControl/paco2n", 40.);
    commonControlParams[1] = ifileC("cardiorespiratoryControl/caco2n", 0.36);

    // [brainAutoParams]
    brainAutoParams[0] = ifile("brainAutoParams/gbpn",0.15);
    brainAutoParams[1] = ifile("brainAutoParams/gbo2",10.);
    brainAutoParams[2] = ifile("brainAutoParams/cvbo2n",0.14);
    brainAutoParams[3] = ifile("brainAutoParams/Ab",20.9);
    brainAutoParams[4] = ifile("brainAutoParams/Bb",92.8);
    brainAutoParams[5] = ifile("brainAutoParams/Cb",10570.);
    brainAutoParams[6] = ifile("brainAutoParams/Db",-5.251);

    // [coroAutoParams]
    coromuscleAutoParams[0] = ifile("coroAutoParams/gho2",35.);
    coromuscleAutoParams[2] = ifile("coroAutoParams/cvho2n",0.11);
    coromuscleAutoParams[4] = ifile("coroAutoParams/khco2",11.11);

    // [muscleAutoParams]
    coromuscleAutoParams[1] = ifile("muscleAutoParams/gmo2",30.);
    coromuscleAutoParams[3] = ifile("muscleAutoParams/cvmo2n",0.155);
    coromuscleAutoParams[5] = ifile("muscleAutoParams/kmco2",142.8);

    // [commonAutoParams]
    commonAutoParams[0] = ifile("commonAutoParams/tauo2A",10.);
    commonAutoParams[1] = ifile("commonAutoParams/tauco2A",20.);

    // [baroAfferentParam]
    baroAfferentParams[0] = ifile("baroAfferentParam/tauzB",6.37);
    baroAfferentParams[1] = ifile("baroAfferentParam/taupB",2.076);
    baroAfferentParams[2] = ifile("baroAfferentParam/fabmin",2.52);
    baroAfferentParams[3] = ifile("baroAfferentParam/fabmax",47.78);
    baroAfferentParams[4] = ifile("baroAfferentParam/pn",92.);
    baroAfferentParams[5] = ifile("baroAfferentParam/kab",11.76);

    // [chemoAfferentParam]
    chemoAfferentParams[0] = ifile("chemoAfferentParam/Ac",600.);
    chemoAfferentParams[1] = ifile("chemoAfferentParam/Bc",10.18);
    chemoAfferentParams[2] = ifile("chemoAfferentParam/kco2",1.);
    chemoAfferentParams[3] = ifile("chemoAfferentParam/ko2",200.);
    chemoAfferentParams[4] = ifile("chemoAfferentParam/kstat",20.);
    chemoAfferentParams[5] = ifile("chemoAfferentParam/kdyn",45.);
    chemoAfferentParams[6] = ifile("chemoAfferentParam/tauCAP",3.5);
    chemoAfferentParams[7] = ifile("chemoAfferentParam/tauccCO2dyn",600.);

    // [lsrAfferentParam]
    lsrAfferentParams[0] = ifile("lsrAfferentParam/gap",12.);
    lsrAfferentParams[1] = ifile("lsrAfferentParam/tauLSR",2.);

    // [efferentSympParams]
    efferentSympParams[0] = ifile("efferentSympParams/fesinf",2.1);
    efferentSympParams[1] = ifile("efferentSympParams/fes0",16.11);
    efferentSympParams[2] = ifile("efferentSympParams/fesmax",60.);
    efferentSympParams[3] = ifile("efferentSympParams/kes",0.0675);
    efferentSympParams[4] = ifile("efferentSympParams/Wbsh",-1.75);
    efferentSympParams[5] = ifile("efferentSympParams/Wbsp",-1.1375);
    efferentSympParams[6] = ifile("efferentSympParams/Wbsv",-1.1375);
    efferentSympParams[7] = ifile("efferentSympParams/Wcsh",1.);
    efferentSympParams[8] = ifile("efferentSympParams/Wcsp",1.56);
    efferentSympParams[9] = ifile("efferentSympParams/Wcsv",1.56);
    efferentSympParams[10] = ifile("efferentSympParams/Wpsh",0.);
    efferentSympParams[11] = ifile("efferentSympParams/Wpsp",-0.4250);
    efferentSympParams[12] = ifile("efferentSympParams/Wpsv",-0.4250);

    // [efferentVagalParams]
    efferentVagalParams[0] = ifile("efferentVagalParams/fevinf",6.3);
    efferentVagalParams[1] = ifile("efferentVagalParams/fev0",3.2);
    efferentVagalParams[2] = ifile("efferentVagalParams/fab0",25.);
    efferentVagalParams[3] = ifile("efferentVagalParams/kev",7.06);
    efferentVagalParams[4] = ifile("efferentVagalParams/Wcv",0.208);
    efferentVagalParams[5] = ifile("efferentVagalParams/Wpv",-0.103);
    efferentVagalParams[6] = ifile("efferentVagalParams/thetav",-0.68);

    // [CNShypoxicResponseParam]
    CNSResponseParams[0] = ifile("CNShypoxixResponseParam/chish",53.);
    CNSResponseParams[1] = ifile("CNShypoxixResponseParam/chisp",6.);
    CNSResponseParams[2] = ifile("CNShypoxixResponseParam/chisv",6.);
    CNSResponseParams[3] = ifile("CNShypoxixResponseParam/Ptildeo2sh",45.);
    CNSResponseParams[4] = ifile("CNShypoxixResponseParam/Ptildeo2sp",30.);
    CNSResponseParams[5] = ifile("CNShypoxixResponseParam/Ptildeo2sv",30.);
    CNSResponseParams[6] = ifile("CNShypoxixResponseParam/kiscsh",6.);
    CNSResponseParams[7] = ifile("CNShypoxixResponseParam/kiscsp",2.);
    CNSResponseParams[8] = ifile("CNShypoxixResponseParam/kiscsv",2.);
    CNSResponseParams[9] = ifile("CNShypoxixResponseParam/tauisc",30.);
    CNSResponseParams[10] = ifile("CNShypoxixResponseParam/taucc",20.);
    CNSResponseParams[11] = ifile("CNShypoxixResponseParam/thetashn",3.6);
    CNSResponseParams[12] = ifile("CNShypoxixResponseParam/thetaspn",13.32);
    CNSResponseParams[13] = ifile("CNShypoxixResponseParam/thetasvn",13.32);
    CNSResponseParams[14] = ifile("CNShypoxixResponseParam/gccsh",1.);
    CNSResponseParams[15] = ifile("CNShypoxixResponseParam/gccsp",1.5);
    CNSResponseParams[16] = ifile("CNShypoxixResponseParam/gccsv",0.);

    // [resistanceControlParams]
    resControlParams[0] = ifile("resistanceControlParams/Grmp",2.47);
    resControlParams[1] = ifile("resistanceControlParams/Grsp",0.695);
    resControlParams[2] = ifile("resistanceControlParams/Grep",1.94);
    resControlParams[3] = ifile("resistanceControlParams/Drmp",2.);
    resControlParams[4] = ifile("resistanceControlParams/Drsp",2.);
    resControlParams[5] = ifile("resistanceControlParams/Drep",2.);
    resControlParams[6] = ifile("resistanceControlParams/fesmin",2.66);
    resControlParams[7] = ifile("resistanceControlParams/taurmp",6.);
    resControlParams[8] = ifile("resistanceControlParams/taursp",6.);
    resControlParams[9] = ifile("resistanceControlParams/taurep",6.);

    // [vuControlParams]
    vuControlParams[0] = ifile("vuControlParams/Gvumv",-58.29);
    vuControlParams[1] = ifile("vuControlParams/Gvusv",-265.4);
    vuControlParams[2] = ifile("vuControlParams/Gvuev",-74.21);
    vuControlParams[3] = ifile("vuControlParams/Dvumv",5.);
    vuControlParams[4] = ifile("vuControlParams/Dvusv",5.);
    vuControlParams[5] = ifile("vuControlParams/Dvuev",5.);
    vuControlParams[6] = ifile("vuControlParams/tauvumv",20.);
    vuControlParams[7] = ifile("vuControlParams/tauvusv",20.);
    vuControlParams[8] = ifile("vuControlParams/tauvuev",20.);

    // [heartControlParams]
    heartControlParams[0] = ifile("heartControlParams/Gemaxlv",0.475);
    heartControlParams[1] = ifile("heartControlParams/Gemaxrv",0.282);
    heartControlParams[2] = ifile("heartControlParams/Demaxlv",2.);
    heartControlParams[3] = ifile("heartControlParams/Demaxrv",2.);
    heartControlParams[4] = ifile("heartControlParams/tauemaxlv",8.);
    heartControlParams[5] = ifile("heartControlParams/tauemaxrv",8.);
    heartControlParams[6] = ifile("heartControlParams/GTs",-0.13);
    heartControlParams[7] = ifile("heartControlParams/DTs",2.);
    heartControlParams[8] = ifile("heartControlParams/tauTs",2.);
    heartControlParams[9] = ifile("heartControlParams/GTv",-0.09);
    heartControlParams[10] = ifile("heartControlParams/DTv",0.2);
    heartControlParams[11] = ifile("heartControlParams/tauTv",1.5);

    if (verbose)
		printCardioControlParameters();
}

void cardiovascularControl::InitialCond(string _file){
	if (verbose)
		cout << " -- cardiovascularControl::saveInitialCond() reading from "
				<< _file << endl;

	GetPot ifile(_file.c_str());
    // Initial Conditions
    stateVarAuto[0] = ifile("ics/xbo2" , 0.);
    stateVarAuto[1] = ifile("ics/xbco2" , 0.);
    stateVarAuto[2] = ifile("ics/xho2" , 0.);
    stateVarAuto[3] = ifile("ics/xhco2" , 0.);
    stateVarAuto[4] = ifile("ics/xmo2" , 0.);
    stateVarAuto[5] = ifile("ics/xmco2" , 0.);
    stateVarBaro = ifile("ics/ptilde", 306.8401);
    stateVarChemo[0] = ifile("ics/phiCO2dyn", 82.2871);
    stateVarChemo[1] = ifile("ics/phiapc", 0.);
    stateVarLSR = ifile("ics/fap",0.);
    stateVarCNS[0] = ifile("ics/deltathetao2sh",0.);
    stateVarCNS[1] = ifile("ics/deltathetaco2sh",0.);
    stateVarCNS[2] = ifile("ics/deltathetao2sp",0.);
    stateVarCNS[3] = ifile("ics/deltathetaco2sp",0.);
    stateVarCNS[4] = ifile("ics/deltathetao2sv",0.);
    stateVarCNS[5] = ifile("ics/deltathetaco2sv",0.);
    stateVarSVResponse[0] = ifile("ics/deltarmp",0.);
    stateVarSVResponse[1] = ifile("ics/deltarsp",0.);
    stateVarSVResponse[2] = ifile("ics/deltarep",0.);
    stateVarSVResponse[3] = ifile("ics/deltavumv",0.);
    stateVarSVResponse[4] = ifile("ics/deltavusv",0.);
    stateVarSVResponse[5] = ifile("ics/deltavuev",0.);
    stateVarSVResponse[6] = ifile("ics/deltaemaxlv",0.);
    stateVarSVResponse[7] = ifile("ics/deltaemaxrv",0.);
    stateVarSVResponse[8] = ifile("ics/deltats",0.);
    stateVarSVResponse[9] = ifile("ics/deltatv",0.);

    setStateVarsVector();
    getAlgebraicRelations();
}

void cardiovascularControl::printCardioControlParameters(){
    cout << "[commonParams]"<<endl;
    for (int i = 0; i<11; i++){
        cout << i << " value: " << commonParams[i] << endl;
    }    

    cout << "[commonControlParams]"<<endl;
    for (int i = 0; i<2; i++){
        cout << i << " value: " << commonControlParams[i] << endl;
    }    

    cout << "[brainAutoParams]"<<endl;
    for (int i = 0; i<7; i++){
        cout << i << " value: " << brainAutoParams[i] << endl;
    }    

    cout << "[coromuscleAutoParams]"<<endl;
    for (int i = 0; i<6; i++){
        cout << i << " value: " << coromuscleAutoParams[i] << endl;
    }   

    cout << "[commonAutoParams]"<<endl;
    for (int i = 0; i<2; i++){
        cout << i << " value: " << commonAutoParams[i] << endl;
    }

    cout << "[varoAfferentParams]"<<endl;
    for (int i = 0; i<6; i++){
        cout << i << " value: " << baroAfferentParams[i] << endl;
    }       

    cout << "[baroAfferentParams]"<<endl;
    for (int i = 0; i<6; i++){
        cout << i << " value: " << baroAfferentParams[i] << endl;
    }  

    cout << "[chemoAfferentParams]"<<endl;
    for (int i = 0; i<8; i++){
        cout << i << " value: " << chemoAfferentParams[i] << endl;
    }  

    cout << "[lsrAfferentParams]"<<endl;
    for (int i = 0; i<2; i++){
        cout << i << " value: " << lsrAfferentParams[i] << endl;
    }

    cout << "[efferentSympParams]"<<endl;
    for (int i = 0; i<13; i++){
        cout << i << " value: " << efferentSympParams[i] << endl;
    }  

    cout << "[efferentVagalParams]"<<endl;
    for (int i = 0; i<7; i++){
        cout << i << " value: " << efferentVagalParams[i] << endl;
    }  

    cout << "[CNShypoxicResponseParam]"<<endl;
    for (int i = 0; i<17; i++){
        cout << i << " value: " << CNSResponseParams[i] << endl;
    }  

    cout << "[resControlParams]"<<endl;
    for (int i = 0; i<10; i++){
        cout << i << " value: " << resControlParams[i] << endl;
    }  

    cout << "[vuControlParams]"<<endl;
    for (int i = 0; i<9; i++){
        cout << i << " value: " << vuControlParams[i] << endl;
    }  

    cout << "[heartControlParams]"<<endl;
    for (int i = 0; i<12; i++){
        cout << i << " value: " << heartControlParams[i] << endl;
    }  

}



void cardiovascularControl::getAlgebraicRelationsAuto(){
    // [AUTOREGULATION]
    double Paco2Diff = (varsTransport[0] - commonControlParams[0] );    
    // equivalent to Ursino's model (change of variable)
    algVarAuto[0] = (1. - brainAutoParams[1] * stateVarAuto[0] + stateVarAuto[1]);
    algVarAuto[1] = (brainAutoParams[3] + (brainAutoParams[4]) / (1. + brainAutoParams[5] * exp(brainAutoParams[6]*log(varsTransport[0]))) ) 
                        / (brainAutoParams[3] + (brainAutoParams[4]) / (1. + brainAutoParams[5] * exp(brainAutoParams[6]*log(commonControlParams[0])))) - 1.;
        

    for (int i=0; i<2;i++){
        algVarAuto[i+2] = (1. - exp( Paco2Diff / coromuscleAutoParams[i+4] ) )
                         / (1. + exp(Paco2Diff / coromuscleAutoParams[i+4] ) );
        algVarAuto[i+5] = (1. + stateVarAuto[2*i + 3] ) / (1. - coromuscleAutoParams[i] * stateVarAuto[2*i + 2]);
    }

    algVarAuto[4] = commonParams[8] / algVarAuto[0];
    algVarAuto[5] *= commonParams[9]; 
    algVarAuto[6] *= commonParams[0];
}


void cardiovascularControl::getAlgebraicRelationsBaro(){
    // [BAROREFLEX]
    algVarBaro = ( baroAfferentParams[2] + baroAfferentParams[3] * exp((stateVarBaro - baroAfferentParams[4]) / baroAfferentParams[5] )) 
            / (1 + exp((stateVarBaro - baroAfferentParams[4]) / baroAfferentParams[5]) ); 

}


void cardiovascularControl::getAlgebraicRelationsChemo(){
    // [CHEMOREFLEX]
    algVarChemo[0] = chemoAfferentParams[0] * (1. - varsTransport[2]) + chemoAfferentParams[1];
    algVarChemo[1] = chemoAfferentParams[2] * (varsTransport[3] - commonControlParams[1]) 
                        * chemoAfferentParams[3] * (1. - exp(-algVarChemo[0] / chemoAfferentParams[3]));
    
    if (varsTransport[3] > commonControlParams[1]){
        algVarChemo[2] = chemoAfferentParams[4] * (1. - exp(-algVarChemo[1] / chemoAfferentParams[4]));
    }else{
        algVarChemo[2] = 0.;
    }
    
    algVarChemo[3] = chemoAfferentParams[5] * (1. - exp(- stateVarChemo[0] / chemoAfferentParams[5] ));
    algVarChemo[4] = max(0., stateVarChemo[1]);
}

void cardiovascularControl::getAlgebraicRelationsLSR(){
    algVarLSR = lsrAfferentParams[0] * Vt;
}

void cardiovascularControl::getAlgebraicRelationsCNS(){
    // [CNS ischemic response]
    for (int i = 0; i<3; i++){
        algVarCNS[i] = CNSResponseParams[i] / (1 + exp((varsTransport[1] - CNSResponseParams[i+3]) / (CNSResponseParams[i+6])));
        algVarCNS[i+3] = CNSResponseParams[i+11] - stateVarCNS[i*2] - CNSResponseParams[i+14] * stateVarCNS[i*2+1];
    }
}

void cardiovascularControl::getAlgebraicRelationsFiring(){
    // [FIRING RATES (sympathetic and vagal)]
    for (int i = 0; i<3; i++){
        algVarFiring[i] = efferentSympParams[0] + (efferentSympParams[1] - efferentSympParams[0]) * exp(efferentSympParams[3] * 
                        (efferentSympParams[i+4] * algVarBaro + efferentSympParams[i+7] * algVarChemo[4] + efferentSympParams[i+10] * stateVarLSR - algVarCNS[i+3]) );

        if (algVarFiring[i] > efferentSympParams[2]){
            algVarFiring[i] = efferentSympParams[2];
        }
        else if (algVarFiring[i] < resControlParams[6]){
            algVarFiring[i] = resControlParams[6];
        }
    }

    algVarFiring[3] = (efferentVagalParams[1] + efferentVagalParams[0] * exp((algVarBaro - efferentVagalParams[2]) / efferentVagalParams[3])) / (1. + exp((algVarBaro - efferentVagalParams[2]) / efferentVagalParams[3])) 
                        + efferentVagalParams[4] * algVarChemo[4] + efferentVagalParams[5] * stateVarLSR - efferentVagalParams[6];

}


void cardiovascularControl::getAlgebraicRelationsSVResponse(){
    // Is it correct? I cannot figure where the python implementation comes from
    for (int i = 0; i < 3; i++){
        algVarSVResponse[i] = resControlParams[i] * log(max(delayedVars[i] - resControlParams[6] + 1., 1.));
        algVarSVResponse[i+3] = resControlParams[i] * log(max(delayedVars[i+3] - resControlParams[6] + 1., 1.));
    }
    
    algVarSVResponse[6] = heartControlParams[0] * log(max(delayedVars[6] - resControlParams[6] + 1. , 1.));
    algVarSVResponse[7] = heartControlParams[1] * log(max(delayedVars[7] - resControlParams[6] + 1. , 1.));
    
    algVarSVResponse[8] = heartControlParams[6] * log(max(delayedVars[8] - resControlParams[6] + 1. , 1.));
    algVarSVResponse[9] = heartControlParams[9] * delayedVars[9];
}


void cardiovascularControl::getAlgebraicRelationscontrolledVars(){
    controlledVars[0] = stateVarSVResponse[0]*(1 + stateVarAuto[5] ) / (1. - coromuscleAutoParams[1] * stateVarAuto[4]) + algVarAuto[6]; // autoregulation has and impact as well
    
    for (int i = 1; i<8; i++){
        controlledVars[i] = stateVarSVResponse[i] + commonParams[i];
    }

    controlledVars[8] = stateVarSVResponse[8] + stateVarSVResponse[9] + T0; // cardiac cycle
    controlledVars[9] = algVarAuto[4]; // affected only by autoregulation   
    controlledVars[10] = algVarAuto[5]; // affected only by autoregulation
}


void cardiovascularControl::getTimeDerivativeAuto(){
    for (int i=0; i<3; i++){
        dstateVarAuto[i*2+1] = 1/commonAutoParams[1] * (-stateVarAuto[i*2+1] + algVarAuto[i+1]);
    }


    dstateVarAuto[0] = 1/commonAutoParams[0] * (-stateVarAuto[0] + (varsTransport[4] - brainAutoParams[2]));
    dstateVarAuto[2] = 1/commonAutoParams[0] * (-stateVarAuto[2] + (varsTransport[5] - coromuscleAutoParams[2]));
    dstateVarAuto[4] = 1/commonAutoParams[0] * (-stateVarAuto[4] + (varsTransport[6] - coromuscleAutoParams[3]));
}


void cardiovascularControl::getTimeDerivativeBaro(){
    dstateVarBaro = 1/baroAfferentParams[1] * (Pb + baroAfferentParams[0] * dPbdt - stateVarBaro);
}


void cardiovascularControl::getTimeDerivativeChemo(){
    dstateVarChemo[0] = 1/chemoAfferentParams[6] * (chemoAfferentParams[7] * varsTransport[7] - stateVarChemo[0]);
    dstateVarChemo[1] = 1/chemoAfferentParams[6] * (algVarChemo[2] + algVarChemo[3] - stateVarChemo[1]);
}

void cardiovascularControl::getTimeDerivativeLSR(){
    dstateVarLSR = 1/lsrAfferentParams[1] * (-stateVarLSR + algVarLSR);
}
    

void cardiovascularControl::getTimeDerivativeCNS(){    
    for (int i = 0; i<3; i++){
        dstateVarCNS[i*2] = 1/CNSResponseParams[9] * (-stateVarCNS[i*2] + algVarCNS[i]);
        dstateVarCNS[i*2 + 1] = 1/CNSResponseParams[10] * (-stateVarCNS[i*2 + 1] + (varsTransport[0]-commonControlParams[0]));
    }
}


void cardiovascularControl::getTimeDerivativeSVResponse(){  
    for (int i = 0; i<3; i++){
        dstateVarSVResponse[i] = 1/resControlParams[i+7] * (-stateVarSVResponse[i] + algVarSVResponse[i]);
        dstateVarSVResponse[i+3] = 1/vuControlParams[i+6] * (-stateVarSVResponse[i+3] + algVarSVResponse[i+3]);
    }

    dstateVarSVResponse[6] = 1/heartControlParams[4] * (-stateVarSVResponse[6] + algVarSVResponse[6]);
    dstateVarSVResponse[7] = 1/heartControlParams[5] * (-stateVarSVResponse[7] + algVarSVResponse[7]);
    dstateVarSVResponse[8] = 1/heartControlParams[8] * (-stateVarSVResponse[8] + algVarSVResponse[8]);
    dstateVarSVResponse[9] = 1/heartControlParams[11] * (-stateVarSVResponse[9] + algVarSVResponse[9]);
}


// overall functions

void cardiovascularControl::getAlgebraicRelations(){
    
    getAlgebraicRelationsAuto(); 
    getAlgebraicRelationsBaro();
    getAlgebraicRelationsChemo();
    getAlgebraicRelationsCNS();
    getAlgebraicRelationsLSR();
    getAlgebraicRelationsFiring();
    getAlgebraicRelationsSVResponse();
    getAlgebraicRelationscontrolledVars();


    //get AlgVars
    algVars.resize(2+algVarAuto.size()+algVarChemo.size()+algVarCNS.size()+algVarFiring.size()+algVarSVResponse.size()+controlledVars.size());
    int count = 0;
    for (int i = 0; i < algVarAuto.size(); i++) {
        algVars[i] = algVarAuto[i]; // global variables
        count++;
    }
    algVars[count] = algVarBaro;
    count++;
    for (int i = 0; i < algVarChemo.size(); i++) {
        algVars[count] = algVarChemo[i];
        count++;
    }
    for (int i = 0; i < algVarCNS.size(); i++) {
        algVars[count] = algVarCNS[i];
        count++;
    }
    algVars[count] = algVarLSR;
    count++;
    for (int i = 0; i < algVarFiring.size(); i++) {
            algVars[count] = algVarFiring[i];
        count++;
    }
    for (int i = 0; i < algVarSVResponse.size(); i++) {
            algVars[count] = algVarSVResponse[i];
        count++;
    }
    for (int i = 0; i < controlledVars.size(); i++) {
            algVars[count] = controlledVars[i];
        count++;
    }

}


void cardiovascularControl::getTimeDerivative(){
        // get algebric relations
        getAlgebraicRelations();
        // get time derivatives
        getTimeDerivativeAuto();
        getTimeDerivativeBaro();
        getTimeDerivativeChemo();
        getTimeDerivativeLSR();
        getTimeDerivativeCNS();
        getTimeDerivativeSVResponse();

        //get dvdt
        dvdt.resize(2+dstateVarAuto.size()+dstateVarChemo.size()+dstateVarCNS.size()+dstateVarSVResponse.size());
        int count = 0;
        for (int i = 0; i < dstateVarAuto.size(); i++) {
            dvdt[i] = dstateVarAuto[i]; // global variables
            count++;
        }
        dvdt[count] = dstateVarBaro;
        count++;
        for (int i = 0; i < dstateVarChemo.size(); i++) {
            dvdt[count] = dstateVarChemo[i];
            count++;
        }
        dvdt[count] = dstateVarLSR;
        count++;
        for (int i = 0; i < dstateVarCNS.size(); i++) {
            dvdt[count]= dstateVarCNS[i];
            count++;
        }
        for (int i = 0; i < dstateVarSVResponse.size(); i++) {
            dvdt[count] = dstateVarSVResponse[i];
            count++;
        }
}


void cardiovascularControl::setStateVars(){

    int count = 0;
    for (int i = 0; i < stateVarAuto.size(); i++) {
        stateVarAuto[i] = stateVars[i]; 
        count++;
    }
    stateVarBaro = stateVars[count];
    count++;
    for (int i = 0; i < stateVarChemo.size(); i++) {
         stateVarChemo[i] = stateVars[count];
        count++;
    }
    stateVarLSR = stateVars[count];
    count++;
    for (int i = 0; i < stateVarCNS.size(); i++) {
        stateVarCNS[i] = stateVars[count];
        count++;
    }
    for (int i = 0; i < stateVarSVResponse.size(); i++) {
        stateVarSVResponse[i] = stateVars[count];
        count++;
    }
}

void cardiovascularControl::setStateVarsVector(){

    stateVars.resize(2+stateVarAuto.size()+stateVarChemo.size()+stateVarCNS.size()+stateVarSVResponse.size());
    int count = 0;
    for (int i = 0; i < stateVarAuto.size(); i++) {
        stateVars[i] = stateVarAuto[i]; 
        count++;
    }
    stateVars[count] = stateVarBaro;
    count++;
    for (int i = 0; i < stateVarChemo.size(); i++) {
        stateVars[count] = stateVarChemo[i];
        count++;
    }
    stateVars[count] = stateVarLSR;
    count++;
    for (int i = 0; i < stateVarCNS.size(); i++) {
        stateVars[count] = stateVarCNS[i];
        count++;
    }
    for (int i = 0; i < stateVarSVResponse.size(); i++) {
        stateVars[count] = stateVarSVResponse[i];
        count++;
    }
}

void cardiovascularControl::output(){   
	
    sampleAlgVars  << scientific << setprecision(12) << time << " ";
    for (int i = 0; i<algVars.size(); i++){
        sampleAlgVars << algVars[i] << " " ;
    }                
    sampleAlgVars << "\n";
    sampleAlgVars.flush();

    setStateVarsVector();
    sampleDiffVars  << scientific << setprecision(12) << time << " ";
    for (int i = 0; i<stateVars.size(); i++){
        sampleDiffVars << stateVars[i] << " " ;
    }
    sampleDiffVars << "\n";
    sampleDiffVars.flush();
}
