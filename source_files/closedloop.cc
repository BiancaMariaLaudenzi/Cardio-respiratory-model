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
#include <sys/types.h>
#include <sys/stat.h>
//#include <experimental/filesystem>
#include <cstdlib>
//#include <direct.h> 
#include <cmath>
#include <functional>
#include <algorithm>
#include "closedloop.h"

using namespace std;
using std::scientific;

void closedloop::init(string ifile, int verbose) {
<<<<<<< HEAD
=======

>>>>>>> origin/main
	// verbose of each system
	CVS.verbose = 0;
	lungs.verbose = 0;
	gas.verbose = 0;
	controlCardio.verbose = 0;
	controlResp.verbose = 0;

	// read parameter file for each system
	if (verbose)
		cout << "Reading closedloop params ..." << endl;
	if (ifile != "NOFILE") {
			readParams(ifile, verbose);
	}
	if (verbose){
		// numerical method
		cout << "method :" << method << endl;
		// [conversions]
		cout << "mmHgDyncm2 :" << mmHgDyncm2 << endl;
		cout << "cmH2OmmHg :" << cmH2OmmHg << endl;
		cout << "LmL :" << LmL << endl;
		// [timeConsts]
		cout << "tIni :" << tIni << endl;
		cout << "tIniGas :" << tIniGas << endl;
		cout << "tIniControl :" << tIniControl << endl;
		cout << "timeStop :" << timeStop << endl;
		cout << "tSampleIni :" << tSampleIni << endl;
		cout << "tSampleEnd :" << tSampleEnd << endl;
		cout << "dtSample :" << dtSample << endl;
		cout << "dt :" << dt << endl;
		// [systems]
		cout << "doLung :" << doLung << endl;
		cout << "doGas :" << doGas << endl;
		cout << "doControlCardio :" << doControlCardio << endl;
		cout << "doControlResp :" << doControlResp << endl;
		// [toSample]
		cout << "roundVal :" << roundVal << endl;
		cout << "tol :" << tol << endl;
		cout << "Done reading closedloop params. " << endl;
	}
		
<<<<<<< HEAD
=======

>>>>>>> origin/main
	// [initialize time]
	time = tIni;
	CVS.time = tIni;
    if(doLung){
		lungs.time = tIni;
		if(doGas){
			gas.time = tIni;
			if(doControlCardio){
				controlCardio.time = tIni;
				if (doControlResp){
					controlResp.time = tIni;
				}
			}
		}
	}

	// initialise each system
	CVS.init(ifileCVS, ifileCommon, pathResults); 
	// assign CVS initial conditions for common parameters
	CVS.params[5] *= cmH2OmmHg;
	// delayed quantities
	timeList.resize(0); 
	timeList.push_back(time);
	list_time = 20.;
	// assign variables from other systems
	if(doLung){
		lungs.init(ifileLungs, ifileCommon, pathResults); 
		if(doGas){
			gas.init(ifileGas, ifileCommon, pathResults);
			// delayed quantities
			delayedGasList_0.resize(0); 
			delayedGasList_1.resize(0); 
			delayedGasList_2.resize(0); 
			delayedGasList_3.resize(0); 
			delayedGasList_0.push_back(gas.AlveoliLungsAlg[5]); 
			delayedGasList_1.push_back(gas.AlveoliLungsAlg[6]); 
			delayedGasList_2.push_back(max(gas.stateVars[26],0.)); 
			delayedGasList_3.push_back(max(gas.stateVars[27],0.)); 
			// set gasTransport variables from other systems
			gas.CommonVars[0]=lungs.vd*LmL;
			gas.CommonVars[1]=lungs.Vdot*LmL; //
			gas.CommonVars[2]=lungs.VAdot*LmL; //
			gas.CommonVars[3]=lungs.vol[3]*LmL;
			for (int i=1; i < 6; i++){
				gas.CommonVars[i+3]=CVS.stateVars[i];
			}
			gas.CommonVars[9]=CVS.sysCircAlg[27];
			for (int i=9; i < 12; i++){
				gas.CommonVars[i+1]=CVS.stateVars[i];
			}
			gas.CommonVars[12]+= CVS.TveinParam[4];
			gas.CommonVars[13]=CVS.stateVars[13];
			for (int i=7; i < 12; i++){
				gas.CommonVars[i+7]=CVS.sysCircAlg[i];
			}
			for (int i=20; i < 25; i++){
				gas.CommonVars[i-1]=CVS.sysCircAlg[i];
			}
			for (int i=26; i < 33; i++){
				gas.CommonVars[i-2]=CVS.sysCircAlg[i];
			}
			for (int i=2; i < 4; i++){
				gas.CommonVars[i+29]=CVS.pulCircAlg[i];
			}
			if (doControlCardio){
				controlCardio.init(ifileControlCardio, ifileCommon, pathResults);
				// delayed quantities
				delayedControlCardioList_0.resize(0); 
				delayedControlCardioList_1.resize(0);
				delayedControlCardioList_2.resize(0);
				delayedControlCardioList_4.resize(0);
				derivativeControlCardioList_0.resize(0);
				derivativeControlCardioList_0.resize(0);
				delayedControlCardioList_0.push_back(controlCardio.algVarFiring[1]);
				delayedControlCardioList_1.push_back(controlCardio.algVarFiring[2]);
				delayedControlCardioList_2.push_back(controlCardio.algVarFiring[0]);
				delayedControlCardioList_4.push_back(controlCardio.algVarFiring[3]);
				derivativeControlCardioList_0.push_back(gas.AlveoliLungsAlg[6]-controlCardio.commonControlParams[0]);
				derivativeControlCardioList_1.push_back(CVS.sysCircAlg[0]);
				// from gas transport
				controlCardio.varsTransport[0] = gas.AlveoliLungsAlg[8];
				controlCardio.varsTransport[1] = gas.AlveoliLungsAlg[7];
				gas.additionalRelations();
				controlCardio.varsTransport[2] = gas.sao2;
				controlCardio.varsTransport[3] = gas.AlveoliLungsAlg[6];
				controlCardio.varsTransport[4] = gas.stateVars[14];
				controlCardio.varsTransport[5] = gas.stateVars[12];
				controlCardio.varsTransport[6] = gas.stateVars[10];
				controlCardio.varsTransport[7] = 0.;
				
				// tidal volume 
				lungVolList.resize(0);  
				timeVolList.resize(0);
				lungVolList.push_back(lungs.v);
				timeVolList.push_back(time);
				Vtidal = 0.;
				controlCardio.Vt = Vtidal;

				// from cardiovascular
				controlCardio.Pb = CVS.sysCircAlg[0];
				controlCardio.dPbdt = 0.;

				// common vars with CVS
				for(int i = 0; i < 6;i++){
					controlCardio.controlledVars[i] = CVS.controlledVars[6-i];
				}
				controlCardio.controlledVars[6] = CVS.controlledVars[9];
				controlCardio.controlledVars[7] = CVS.controlledVars[10];
				controlCardio.controlledVars[8] = CVS.controlledVars[0];
				controlCardio.controlledVars[9] = CVS.controlledVars[8];
				controlCardio.controlledVars[10] = CVS.controlledVars[7];
				
				if (doControlResp){
					controlResp.init(ifileControlResp, ifileCommon, pathResults);
					// common vars with CardioControl
					delayedControlRespList_0.resize(0);
					delayedControlRespList_1.resize(0);
					delayedControlRespList_0.push_back(controlCardio.varsTransport[0]);
					delayedControlRespList_1.push_back(controlCardio.algVarChemo[4]);
					controlResp.delayedVars[0] = controlCardio.varsTransport[0];
					controlResp.delayedVars[1] = controlCardio.algVarChemo[4];
					
				}
			}
		}
	}
	// save initial conditions for state variables and algebric variables
	CVS.output();
	if (doLung){
		lungs.output();
		if (doGas){
			gas.output();
			if (doControlCardio){
				controlCardio.output();
				if (doControlResp){
					controlResp.output();
					}
			}
		}
	}
		// to built vector of state variables
	nVars = CVS.stateVars.size();
	if (doLung){
		nVars += lungs.vol.size();
		if (doGas){
			nVars += gas.stateVars.size();
			if (doControlCardio){
				nVars += controlCardio.stateVars.size();
				if (doControlResp){
					nVars += controlResp.stateVars.size();
				}
			}
		}
	}	
	x.resize(nVars);
	// ensure that gas transport starts after delayed quantities timing
	if (doGas){
		if (tIniGas < gas.DelayedGasParam[0] or tIniGas < gas.DelayedGasParam[1]){
			cout << "You can not start transport before delayed quantities timing" << endl;
			cout << "Change tIniTransport in closeedloopmodel input file" <<endl;
		}
	}

}

void closedloop::exchangeInfo(){
	if(doLung){
		CVS.params[5] = lungs.P[4] * cmH2OmmHg;
		if(doGas){
			// set gasTransport variables from other systems
			gas.CommonVars[0]=lungs.vd*LmL;
			gas.CommonVars[1]=lungs.Vdot*LmL; //
			gas.CommonVars[2]=lungs.VAdot*LmL; //
			gas.CommonVars[3]=lungs.vol[3]*LmL;
			for (int i=1; i < 6; i++){
				gas.CommonVars[i+3]=CVS.stateVars[i];
			}
			gas.CommonVars[9]=CVS.sysCircAlg[27];
			for (int i=9; i < 12; i++){
				gas.CommonVars[i+1]=CVS.stateVars[i];
			}
			gas.CommonVars[12]+= CVS.TveinParam[4];
			gas.CommonVars[13]=CVS.stateVars[13];
			for (int i=7; i < 12; i++){
				gas.CommonVars[i+7]=CVS.sysCircAlg[i];
			}
			for (int i=20; i < 25; i++){
				gas.CommonVars[i-1]=CVS.sysCircAlg[i];
			}
			for (int i=26; i < 33; i++){
				gas.CommonVars[i-2]=CVS.sysCircAlg[i];
			}
			for (int i=2; i < 4; i++){
				gas.CommonVars[i+29]=CVS.pulCircAlg[i];
			}
			if (doControlCardio){
				// from gas transport
				controlCardio.varsTransport[0] = gas.AlveoliLungsAlg[8];
				controlCardio.varsTransport[1] = gas.AlveoliLungsAlg[7];
				gas.additionalRelations();
				controlCardio.varsTransport[2] = gas.sao2;
				controlCardio.varsTransport[3] = gas.AlveoliLungsAlg[6];
				controlCardio.varsTransport[4] = gas.stateVars[14];
				controlCardio.varsTransport[5] = gas.stateVars[12];
				controlCardio.varsTransport[6] = gas.stateVars[10];
				
				// tidal volume from lungMechanics
				controlCardio.Vt = Vtidal;

				// from cardiovascular
				controlCardio.Pb = CVS.sysCircAlg[0];

				// change CVS params
				for(int i = 0; i < 6;i++){
					CVS.controlledVars[6-i]=controlCardio.controlledVars[i];
				}
				CVS.controlledVars[9] = controlCardio.controlledVars[6];
				CVS.controlledVars[10] = controlCardio.controlledVars[7];
				CVS.controlledVars[0] = controlCardio.controlledVars[8];
				CVS.controlledVars[8] = controlCardio.controlledVars[9];
				CVS.controlledVars[7] = controlCardio.controlledVars[10];
				if (doControlResp){
					// change lungMechanics params
					lungs.resp[0] = controlResp.algVars[1]; 
					lungs.resp[4] = controlResp.algVars[0];  
				}
			}
		}
	}
}

vector<double> closedloop::getStateTimeDerivative(double time){
    vector<double> dxdt(nVars,0);

	exchangeInfo();

	CVS.getTimeDerivative();
	int i0 = 0;
	int i1 = CVS.stateVars.size();
	for (int i =0; i <i1-i0;i++){
		dxdt[i+i0]=CVS.dvdt[i];
	}
	if (doLung){
		lungs.getTimeDerivative();
		i0 = CVS.stateVars.size();
		i1 = i0 + lungs.vol.size();
		for (int i =0; i <i1-i0;i++){
			dxdt[i+i0]=lungs.dvdt[i];
		}
		if (doGas && time > tIniGas) {
			gas.getTimeDerivative();
			i0 = CVS.stateVars.size() + lungs.vol.size();
			i1 = i0 + gas.stateVars.size();
			for (int i =0; i <i1-i0;i++){
				dxdt[i+i0]=gas.dvdt[i];
			}
			if (doControlCardio && time > tIniControl){
				controlCardio.getTimeDerivative();
				i0 = CVS.stateVars.size() + lungs.vol.size()+ gas.stateVars.size();
				i1 = i0 + controlCardio.stateVars.size();
				for (int i =0; i <i1-i0;i++){
					dxdt[i+i0]=controlCardio.dvdt[i];
				}
				if (doControlResp && time > tIniControl){
					controlResp.getTimeDerivative();
					i0 = CVS.stateVars.size() + lungs.vol.size()+ gas.stateVars.size()+controlCardio.stateVars.size();
					i1 = i0 + controlResp.stateVars.size();
					for (int i =0; i <i1-i0;i++){
						dxdt[i+i0]=controlResp.dvdt[i];
					}
				}
			}
		}
	}
	return dxdt;
}

vector<double> closedloop::setStateFromVariables(){
	vector<double> x0(0,0);

	// delayed quantities
	timeList.push_back(time); // append time	
	if (time > list_time){
		// remove first element
		timeList.pop_front(); 
	}

	int i0 = 0;
	x0.insert(x0.begin()+i0, CVS.stateVars.begin(), CVS.stateVars.end());

	if (doLung){

		i0 = CVS.stateVars.size();
		x0.insert(x0.begin()+i0, lungs.vol.begin(), lungs.vol.end());

		if (doGas){
			// delayed quantities
			delayedGasList_0.push_back(gas.TissueVenousGasAlg[0]); 
			delayedGasList_1.push_back(gas.TissueVenousGasAlg[1]); 
			delayedGasList_2.push_back(max(gas.TissueVenousGasAlg[2],0.)); 
			delayedGasList_3.push_back(max(gas.TissueVenousGasAlg[3],0.)); 
			if (time > list_time){ 
				// remove first element
				delayedGasList_0.pop_front();
				delayedGasList_1.pop_front();
				delayedGasList_2.pop_front();
				delayedGasList_3.pop_front();
			}
			// interpolation for delayed quantities
			if (time > gas.DelayedGasParam[0]){
				gas.TissueVenousGasAlg[0] = interpolate(timeList,delayedGasList_0, time-gas.DelayedGasParam[0], false); 
				gas.TissueVenousGasAlg[1] = interpolate(timeList,delayedGasList_1, time-gas.DelayedGasParam[0], false);
			}
			else{
				gas.TissueVenousGasAlg[0] = gas.AlveoliLungsAlg[5]; 
				gas.TissueVenousGasAlg[1] = gas.AlveoliLungsAlg[6]; 
			}
			if (time > gas.DelayedGasParam[1]){
				gas.TissueVenousGasAlg[2] = interpolate(timeList,delayedGasList_2, time-gas.DelayedGasParam[1], false); // 0
				gas.TissueVenousGasAlg[3] = interpolate(timeList,delayedGasList_3, time-gas.DelayedGasParam[1], false); // 1
			}
			else{
				gas.TissueVenousGasAlg[2] = gas.stateVars[26]; 
				gas.TissueVenousGasAlg[3] = gas.stateVars[27]; 
			}	

			i0 = CVS.stateVars.size() + lungs.vol.size();
			x0.insert(x0.begin()+i0, gas.stateVars.begin(), gas.stateVars.end());

			if (doControlCardio){

				// delayed quantities
				delayedControlCardioList_0.push_back(controlCardio.algVarFiring[1]);
				delayedControlCardioList_1.push_back(controlCardio.algVarFiring[2]); 
				delayedControlCardioList_2.push_back(controlCardio.algVarFiring[0]); 
				delayedControlCardioList_4.push_back(controlCardio.algVarFiring[3]);
				derivativeControlCardioList_0.push_back(gas.AlveoliLungsAlg[6]-controlCardio.commonControlParams[0]);
				derivativeControlCardioList_1.push_back(CVS.sysCircAlg[0]);
				if (time > list_time) {
					delayedControlCardioList_0.pop_front();
					delayedControlCardioList_1.pop_front();
					delayedControlCardioList_2.pop_front();
					delayedControlCardioList_4.pop_front();
					derivativeControlCardioList_0.pop_front();
					derivativeControlCardioList_1.pop_front();
				}
				// interpolation for delayed quantities
				if (time > controlCardio.resControlParams[4]){
					controlCardio.delayedVars[0] = interpolate(timeList,delayedControlCardioList_0, time-controlCardio.resControlParams[4], false);
					controlCardio.delayedVars[1] = interpolate(timeList,delayedControlCardioList_0, time-controlCardio.resControlParams[4], false); 
					controlCardio.delayedVars[2] = interpolate(timeList,delayedControlCardioList_0, time-controlCardio.resControlParams[4], false); 
				}
				else {
					controlCardio.delayedVars[0] = 0.;
					controlCardio.delayedVars[1] = 0.;
					controlCardio.delayedVars[2] = 0.;
				}
				if (time > controlCardio.vuControlParams[4]){
					controlCardio.delayedVars[3] = interpolate(timeList,delayedControlCardioList_1, time-controlCardio.vuControlParams[4], false); 
					controlCardio.delayedVars[4] = interpolate(timeList,delayedControlCardioList_1, time-controlCardio.vuControlParams[4], false); 
					controlCardio.delayedVars[5] = interpolate(timeList,delayedControlCardioList_1, time-controlCardio.vuControlParams[4], false);
				}
				else {
					controlCardio.delayedVars[3] = 0.;
					controlCardio.delayedVars[4] = 0.;
					controlCardio.delayedVars[5] = 0.;
				}
				if (time > controlCardio.heartControlParams[2]){
					controlCardio.delayedVars[6] = interpolate(timeList,delayedControlCardioList_2, time-controlCardio.heartControlParams[2], false);
					controlCardio.delayedVars[7] = interpolate(timeList,delayedControlCardioList_2, time-controlCardio.heartControlParams[2], false);  
				}
				else {
					controlCardio.delayedVars[6] = 0.;
					controlCardio.delayedVars[7] = 0.;
				}
				if (time > controlCardio.heartControlParams[7]){
					controlCardio.delayedVars[8] = interpolate(timeList,delayedControlCardioList_2, time-controlCardio.heartControlParams[7], false);
				}
				else {
					controlCardio.delayedVars[8] = 0.;
				}
				if (time > controlCardio.heartControlParams[10]){
					controlCardio.delayedVars[9] = interpolate(timeList,delayedControlCardioList_4, time-controlCardio.heartControlParams[10], false);
				}
				else {
					controlCardio.delayedVars[9] = 0.;
				}
				if (derivativeControlCardioList_1.size()>=2){
					controlCardio.dPbdt = (derivativeControlCardioList_1[derivativeControlCardioList_1.size() - 1]-derivativeControlCardioList_1[derivativeControlCardioList_1.size() - 2])/dt;
				}
				else{
					controlCardio.dPbdt = 0.;
				}
				if (derivativeControlCardioList_0.size()>=2){
					controlCardio.varsTransport[7] = (derivativeControlCardioList_0[derivativeControlCardioList_0.size() - 1]-derivativeControlCardioList_0[derivativeControlCardioList_0.size() - 2])/dt;
				}
				else{
					controlCardio.varsTransport[7] = 0.;
				}

				// // tidal volume 
				lungVolList.push_back(lungs.v);
				timeVolList.push_back(time);
				vector<double> var1(0);
<<<<<<< HEAD
				for ( int i = 0 ; i<timeVolList.size();i++){
=======
				for ( int i; i<timeVolList.size();i++){
>>>>>>> origin/main
					if (timeVolList[i] >= time-lungs.T){
						var1.push_back(lungVolList[i]);
					}
				}
				Vtidal = max_element(var1.begin(),var1.end()) - min_element(var1.begin(),var1.end());
				lungVolList.pop_front();
				timeVolList.pop_front();



				
				i0 = CVS.stateVars.size() + lungs.vol.size()+gas.stateVars.size();
				x0.insert(x0.begin()+i0, controlCardio.stateVars.begin(), controlCardio.stateVars.end());

				if (doControlResp){

					// delayed quantities
					delayedControlRespList_0.push_back(controlCardio.varsTransport[0]);
					delayedControlRespList_1.push_back(controlCardio.algVarChemo[4]);
					if (time > list_time) {
						delayedControlRespList_0.pop_front();
						delayedControlRespList_1.pop_front();
					}
					// interpolation for delayed quantities
					if (time > controlResp.Params[8]){
						controlResp.delayedVars[0] = interpolate(timeList,delayedControlRespList_0, time-controlResp.Params[8], false);
					}
					else {
						controlResp.delayedVars[0] = 0.;
					}
					if (time > controlResp.Params[0]){
						
						controlResp.delayedVars[1] = interpolate(timeList,delayedControlRespList_1, time-controlResp.Params[0], false);
					}
					else {
						controlResp.delayedVars[1] = 0.;
					}

					i0 = CVS.stateVars.size() + lungs.vol.size()+gas.stateVars.size()+controlCardio.stateVars.size();
					x0.insert(x0.begin()+i0, controlResp.stateVars.begin(), controlResp.stateVars.end());
				}
			}
		}
	}
	return x0;
}

void closedloop::setVariablesFromState(double time, const vector<double>& x0){
	int  i0 = 0;
	int i1 = CVS.stateVars.size();
	CVS.stateVars = vector<double>(x0.begin() + i0, x0.begin() + i1);
	CVS.getAlgebraicRelations();
	if (doLung){
		i0 = CVS.stateVars.size();
		i1 = i0 + lungs.vol.size();		
		lungs.vol = vector<double>(x0.begin() + i0, x0.begin() + i1);
		lungs.getAlgebraicRelations();
		if (doGas){
			i0 = CVS.stateVars.size() + lungs.vol.size();
			i1 = i0 + gas.stateVars.size();
			gas.stateVars = vector<double>(x0.begin() + i0, x0.begin() + i1);
			gas.getAlgebraicRelations();
			if (doControlCardio){
				i0 = CVS.stateVars.size() + lungs.vol.size() + gas.stateVars.size();
				i1 = i0 + controlCardio.stateVars.size();
				controlCardio.stateVars = vector<double>(x0.begin() + i0, x0.begin() + i1);
				controlCardio.setStateVars();
				controlCardio.getAlgebraicRelations();
				if (doControlResp){
					i0 = CVS.stateVars.size() + lungs.vol.size() + gas.stateVars.size() + controlCardio.stateVars.size();
					i1 = i0 + controlResp.stateVars.size();
					controlResp.stateVars = vector<double>(x0.begin() + i0, x0.begin() + i1);
					controlResp.getAlgebraicRelations();
				}
			}
		}
	}
	exchangeInfo();
}

void closedloop::Euler(double dt, double time, const vector<double>& x0, function<vector<double>(double)> fun){
	vector <double> dvdt = fun(time);
	for (int i = 0; i < nVars; i++) {
		x[i] += dt * dvdt[i];
	}
} 

void closedloop::RK4(double dt, double time, const vector<double>& x0, function<vector<double>(double)> fun){
	vector<double> k1 = fun(time);
	vector <double> x1(nVars,0);
	for (int i = 0; i < nVars; i++) {
		x1[i] = x0[i] + 0.5 * dt * k1[i];
	}
	vector<double> k2 = fun(time + 0.5 * dt);
	vector <double> x2(nVars,0);
	for (int i = 0; i < nVars; i++) {
		x2[i] = x1[i] + 0.5 * dt * k2[i];
	}
	vector<double> k3 = fun(time + 0.5 * dt);
	vector <double> x3(nVars,0);
	for (int i = 0; i < nVars; i++) {
		x3[i] = x2[i] + dt * k3[i];
	}
	vector<double> k4 = fun(time + dt);
	for (int i = 0; i < nVars; i++) {
		x[i] = x0[i] + (dt / 6.) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
	}
}

void closedloop::timeLoop(int verbose){
	cout << "Start simulation time:: tIni = " << tIni << " s." << endl;
	for(int iter = 0; iter < 1000000000; iter++){
		iT = iter;
		if(time+dt>timeStop){
			dt = timeStop - time;
		}
		if(abs(time - timeStop) <= 1e-6){
			break;
		}
		// update
		x = setStateFromVariables();
		auto Fun = [this](double t) {
			return this -> getStateTimeDerivative(t);
		};
		if (method == 0) {		
			Euler(dt, time, x, Fun);
		}		
		if (method == 1) {
			RK4(dt, time, x, Fun);
		}
		// update times
		time += dt;                 //global
		CVS.time += dt;
		if (doLung){
			lungs.time += dt;
			if (doGas){
				gas.time += dt;
				if (doControlCardio){
					controlCardio.time += dt;
					if (doControlResp){
						controlResp.time += dt;
					}
				}
			}
		}
		setVariablesFromState(time, x);

		// save each dtSample
		int sampleNow = 0;
		if (dtSample == 0.) {
			sampleNow = 1;
		}
		else {
			if ((time >= tSampleIni - tol) && (time <= tSampleEnd + tol)) {
				double tAux = (round(time / (dtSample))) * dtSample;
				double tDiff = round(abs(tAux - time) * roundVal) / roundVal;
				if (tDiff < tol) {
					sampleNow = 1;
				}
			}
		}
		if (sampleNow == 1) {
			if (verbose){
				cout << "time = "<< time << " s" << endl;
			}

			CVS.output();
			if (doLung){
				lungs.output();
				if (doGas){
					gas.output();
					if (doControlCardio){
						controlCardio.output();
						if (doControlResp){
							controlResp.output();
							}
					}
				}
			}
		}			

		// correct time each dtSample
		time = double(int(round(time/dt)))*dt;
	}
	cout << "Results sampled from " << tSampleIni << " s to " << tSampleEnd << " s, each " << dtSample << " s." << endl;
	cout << "Final simulation time:: tEnd = " << time << " s." << endl;
}


void closedloop::readParams(string _file, int verbose) {
	if (verbose)
		cout << " -- closedloop::readClosedloopParameters() reading from "
				<< _file << endl;

	GetPot ifile(_file.c_str());

	pathResults = ifile("folders/pathResults","./test_results/");

	int status;
	status = mkdir(pathResults.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); //ubuntu
	//status = mkdir(pathResults.c_str()); //windows

	// #ifdef _WIN32
	// 	status = mkdir(pathResults.c_str());  //windows
	// #else
	// 	status = mkdir(pathResults.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);  //ubuntu
	// #endif
	
    if (status != 0 and verbose)
		cout << "Model1d::init() could not create folder for results "
				<< endl;

	// define files for each system 
	ifileCommon = ifile("commonParams/commonParams", "noCommonFile");
	ifileCVS = ifile("cardioParams/cardioParams", "noCVSfile");
	ifileLungs = ifile("lungParams/lungParams", "noLungMechaniscFile");
	ifileGas = ifile("gasParams/gasParams", "noTransportFile");
	ifileControlCardio = ifile("controlCardioParams/controlCardioParams", "noCardioControlFile");
	ifileControlResp = ifile("controlRespParams/controlRespParams", "noRespControlFile");
	// [method]
	method = ifile("methods/method", 1);
	// [conversions]
	mmHgDyncm2 = ifile("conversions/mmHgDyncm2", 0.);
	cmH2OmmHg = ifile("conversions/cmH2OmmHg", 0.);
	LmL = ifile("conversions/LmL", 0.);
	// [timeConsts]
	tIni = ifile("timeParams/tIni", 0.);
	tIniGas = ifile("timeParams/tIniGas", 0.);
	tIniControl = ifile("timeParams/tIniControl", 0.);
	timeStop = ifile("timeParams/timeStop", 0.);
	tSampleIni = ifile("timeParams/tSampleIni", 0.);
	tSampleEnd = ifile("timeParams/tSampleEnd", 0.);
	dtSample = ifile("timeParams/dtSample", 0.);
	dt = ifile("timeParams/dt", 0.);
	// [systems]
	doLung = ifile("systems/doLung", 0);
	doGas = ifile("systems/doGas", 0);
	doControlCardio = ifile("systems/doControlCardio", 0);
	doControlResp = ifile("systems/doControlResp", 0);
	// [toSample]
	roundVal = ifile("toSample/roundVal", 0.);
	tol = ifile("toSample/tol", 0.);
}

double closedloop::interpolate(const deque<double>& xData, const deque<double>& yData, double x0, bool extrapolate) {
    size_t size = xData.size();

	if (size != yData.size()){
        cout << "Interpolation Error: The datasets must have the same size." << endl;
        exit(-1);
	}

    if (size < 2) {
        cout << "Interpolation Error: The dataset must contain at least two points for interpolation." << endl;
        exit(-1);
    }

    auto it = lower_bound(xData.begin(), xData.end(), x0);
    size_t index = distance(xData.begin(), it);

    if (index == 0)
        index = 1; // Handle the case when x0 is smaller than the smallest value in xData

    if (!extrapolate && (x0 < xData.front() || x0 > xData.back())) {
        return yData.back(); // Without extrapolation, return the last known value beyond the boundaries
    }

    double xL = xData[index - 1], yL = yData[index - 1], xR = xData[index], yR = yData[index];

    return yL + ((x0 - xL) / (xR - xL)) * (yR - yL);
<<<<<<< HEAD
}

=======
}
>>>>>>> origin/main
