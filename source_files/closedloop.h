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
#include <deque>
#include <cmath>
#include <functional>
#include "cardiovascular.h"
#include "lungMechanics.h"
#include "transport.h"
#include "cardiovascularControl.h"
#include "respiratoryControl.h"

using namespace std;


/**
 * @brief Closed-loop
*/

class closedloop
{
  public:

    cardiovascular CVS; // declare CVS of type cardiovascular
	  lungMechanics lungs; // declare lungs of type lungMechanics
    gasTransport gas; 
    cardiovascularControl controlCardio;
    respiratoryControl controlResp;

    // general inputs
    int iT;
    int method;

    // systems
    int doLung;
    int doGas;
    int doControlCardio;
    int doControlResp;

    // time constants
    double tIni;
    double tIniGas;
    double tIniControl;
    double timeStop;
    double tSampleIni;
    double tSampleEnd;
    double dtSample;
    double time;
    double dt;
    double list_time;

    // conversions
    double mmHgDyncm2;
    double cmH2OmmHg;
    double LmL;

    // total number of variables
    int nVars;
    // vector of solutions
    vector<double> x;

    //to sample the solution
    double roundVal;
    double tol;

    // input files
    string ifileCVS;
    string ifileCommon;
    string ifileLungs;
    string ifileGas;
    string ifileControlCardio;
    string ifileControlResp;

    // path to save results
    string pathResults;

    // Delayed quantities
    deque<double> timeList;

    deque<double> delayedGasList_0; // cao2
    deque<double> delayedGasList_1; // caco2
    deque<double> delayedGasList_2; // cvo2
    deque<double> delayedGasList_3; // cvco2

    deque<double> delayedControlCardioList_0; // resistance
    deque<double> delayedControlCardioList_1; // unstressed vol
    deque<double> delayedControlCardioList_2; // emax and T sympathetic
    deque<double> delayedControlCardioList_4; // T vagal
    
    deque<double> derivativeControlCardioList_0; // caco2
    deque<double> derivativeControlCardioList_1; // psa

    deque<double> delayedControlRespList_0; // paco2
    deque<double> delayedControlRespList_1; // fapc

    // tidal volume
    deque<double> lungVolList;
    deque<double> timeVolList;

    //tidal volume 
    double Vtidal;

    // methods
    void init(string ifile, int verbose);
    void readParams(string _file,  int verbose);
    vector<double> getStateTimeDerivative(double time);
    vector<double> setStateFromVariables();
    void setVariablesFromState(double time, const vector<double>& x0);
    void Euler(double dt, double time, const vector<double>& x0, function<vector<double>(double)> fun);
    void RK4(double dt, double time, const vector<double>& x0, function<vector<double>(double)> fun);
    void exchangeInfo();
    void timeLoop(int verbose);
    double interpolate(const deque<double>& xData, const deque<double>& yData, double x0, bool extrapolate);
};