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
 * @brief Cardiovascular model
 * 
 * @note - Vascular model modified from the one proposed in:
 *
 * A Albanese, L Cheng, M Ursino & N. W. Chbat.
 * An integrated mathematical model of the human cardiopulmonary system: model development.
 * American Journal of Physiology-Heart and Circulatory Physiology 2016 310, H899-H921.   (A)
 *
 * Modifications -> 
 * - ODE variables: Volumes instead Pressures; only 2-elements Windkessel models
 *
 * - Heart valve model proposed in:
 *
 * JP Mynard, MR Davidson, DJ Penny, and JJ Smolich. 
 * A simple, versatile valve model for use in the lumped parameter and one-dimensional cardiovascular models. 
 * International Journal for Numerical Methods in the Biomedical Engineering, 28(6-7):626–641, 2012.    (M)
 *
 * - Ventricles model proposed in:
 *
 * Mauro Ursino and Elisa Magosso. 
 * Acute cardiovascular response to isocapnic hypoxia. I. A mathematical model. 
 * American Journal of Physiology-Heart and Circulatory Physiology, 279(1):H149-H165, July 2000.   (U)
 *
 * - Atria model proposed in:
 *
 * FY Liang, S Takagi, R Himeno, and H Liu. 
 * Biomechanical characterization of ventricular–arterial coupling during aging: a multi-scale model study. 
 * Journal of biomechanics, 42(6):692–704, 2009.    (L)
 *
 * Units:
 *      # [compliances] : ml / mmHg 
 *      # [resistances] : mmHg s / ml
 *      # [volumes] : ml
 *      # [areas] : cm^2
*/

class cardiovascular
{
  public:
    /// time 
    double time;
    int verbose;

    //-----------------------------------------------------------------------------------
    // State Variables
    //-----------------------------------------------------------------------------------

    /**
     * @brief Vector with 29 components that contains state variables
     * 
     * Systemic arteries circulation
     * double vsa;     [0]  // volume in the systemic arteries
     * 
     * Systemic peripheral circulation
     * double vep;     [1] // volume in the extrasplanchnic peripheral compartment
     * double vsp;     [2] // volume in the splanchnic peripheral compartment
     * double vmp;     [3] // volume in the skeletal muscle peripheral compartment
     * double vhp;     [4] // volume in the coronary peripheral compartment
     * double vbp;     [5] // volume in the brain peripheral compartment
     * 
     * Systemic venous circulation
     * double vev;     [6] // volume in the extrasplanchnic veins
     * double vsv;     [7] // volume in the splanchnic veins
     * double vmv;     [8] // volume in the skeletal muscle veins 
     * double vhv;     [9] // volume in the coronary veins
     * double vbv;     [10] // volume in the brain veins
     * double vtv;     [11] // volume in the thoracic veins
     * 
     * Pulmonary circulation
     * double vpa;     [12] // volume in the pulmonary artery
     * double vpp;     [13] // volume in the pulmonary peripheral circulation
     * double vpv;     [14] // volume in the pulmonary veins
     * double vps;     [15] // volume in the pulmonary shunt
     * 
     * Heart chambers
     * double vla;     [16] // volume in the left atrium
     * double vlv;     [17] // volume in the left ventricle
     * double epsilon; [18] // variable that evolves the fraction of cardiac cycle u(t)
     * double vra;     [19] // volume in the right atrium
     * double vrv;     [20] // volume in the right ventricle
     * 
     * Heart valves
     * double qTV;     [21] // flow through tricuspid valve
     * double xiTV;    [22] // opening state of tricuspid valve
     * double qPV;     [23] // flow through pulmonary valve
     * double xiPV;    [24] // opening state of pulmonary valve
     * double qMV;     [25] // flow through mitral valve
     * double xiMV;    [26] // opening state of mitral valve
     * double qAV;     [27] // flow through aortic valve
     * double xiAV;    [28] // opening state of aortic valve
    */
    vector<double> stateVars = vector<double>(29);

    /**
     * @brief Vector of 29 components that contains the derivate of the state variables with respect to time (dvdt)
     * 
     * Systemic arteries circulation
     * double dvsadt;      [0] // dt volume in the systemic arteries
     * Systemic peripheral circulation
     * double dvepdt;      [1] // dt volume in the extrasplanchnic peripheral compartment
     * double dvspdt;      [2] // dt volume in the splanchnic peripheral compartment
     * double dvmpdt;      [3] // dt volume in the skeletal muscle peripheral compartment
     * double dvhpdt;      [4] // dt volume in the coronary peripheral compartment
     * double dvbpdt;      [5] // dt volume in the brain peripheral compartment
     * Systemic venous circulation
     * double dvevnewdt;      [6] // dt volume in the extrasplanchnic veins
     * double dvsvnewdt;      [7] // dt volume in the splanchnic veins
     * double dvmvnewdt;      [8] // dt volume in the skeletal muscle veins
     * double dvhvdt;      [9] // dt volume in the coronary veins
     * double dvbvdt;      [10] // dt volume in the brain veins
     * double dvtvdt;      [11] // dt volume in the thoracic veins
     * Pulmonary circulation
     * double dvpadt;      [12] // dt volume in the pulmonary artery
     * double dvppdt;      [13] // dt volume in the pulmonary peripheral circulation
     * double dvpvdt;      [14] // dt volume in the pulmonary veins
     * double dvpsdt;      [15] // dt volume in the pulmonary shunt
     * Heart chambers
     * double dvladt;      [16] // dt volume in the left atrium
     * double dvlvdt;      [17] // dt volume in the left ventricle
     * double depsilondt;  [18] // dt variable that evolves the fraction of cardiac cycle u(t)
     * double dvradt;      [19] // dt volume in the right atrium
     * double dvrvdt;      [20] // dt volume in the right ventricle
     * Heart valves
     * double dqTVdt;      [21] // dt flow through tricuspid valve
     * double dxiTVdt;     [22] // dt opening state of tricuspid valve
     * double dqPVdt;      [23] // dt flow through pulmonary valve
     * double dxiPVdt;     [24] // dt opening state of tricuspid valve
     * double dqMVdt;      [25] // dt flow through mitral valve
     * double dxiMVdt;     [26] // dt opening state of miral valve
     * double dqAVdt;      [27] // dt flow through aortic valve
     * double dxiAVdt;     [28] // dt opening state of aortic valve
    */
    vector<double> dvdt = vector<double>(29);
  
    //-----------------------------------------------------------------------------------
    // Parameters
    //-----------------------------------------------------------------------------------

    /**
      * @brief Vector of 14 components that contains systemic and pulmonary compliances
      * 
      * double csa;       [0] // systemic arterial compliance
      * double cep;       [1] // extra-splanchnic peripheral compliance
      * double csp;       [2] // splanchnic peripheral compliance
      * double cmp;       [3] // skeletal muscle peripheral compliance
      * double chp;       [4] // coronary peripheral compliance
      * double cbp;       [5] // brain peripheral compliance
      * double cev;       [6] // extra-splanchnic venous compliance
      * double csv;       [7] // splanchnic venous compliance
      * double cmv;       [8] // skeletal muscle venous compliance
      * double chv;       [9] // coronary venous compliance
      * double cbv;       [10] // brain venous compliance
      * double cpa;       [11] // pulmonary artery compliance
      * double cpp;       [12] // pulmonary peripheral compliance
      * double cpv;       [13] // pulmonary venous compliance
    */
    vector<double> Cparam = vector<double>(14);  

    /**
     * @brief Vector of 14 components that contains systemic and pulmonary unstressed volumes
     * 
     * double vusa;      [0] // systemic arterial unstressed volume
     * double vuep;      [1] // extra-splanchnic peripheral unstressed volume
     * double vusp;      [2] // splanchnic peripheral unstressed volume
     * double vump;      [3] // skeletal muscle peripheral unstressed volume
     * double vuhp;      [4] // coronary peripheral unstressed volume
     * double vubp;      [5] // brain peripheral unstressed volume
     * double vuev0;     [6] // extra-splanchnic venous unstressed volume (nominal)   // controlled in cardiovascularControl 
     * double vusv0;     [7] // splanchnic venous unstressed volume (nominal)         // controlled in cardiovascularControl 
     * double vumv0;     [8] // skeletal muscle venous unstressed volume (nominal)    // controlled in cardiovascularControl 
     * double vuhv;      [9] // coronary venous unstressed volume
     * double vubv;      [10] // brain venous unstressed volume
     * double vupa;      [11] // pulmonary artery unstressed volume
     * double vupp;      [12] // pulmonary peripheral nstressed volume
     * double vupv;      [13] // pulmonary venous unstressed volume
     */
    vector<double> VuParam = vector<double>(14);      

    /**
     *  @brief Vector of 14 components that contains systemic and pulmonary resistances
     * 
     * double Rp;       [0] // resistance from systemic artery to peripheral compartments (fixed at 0.001 mmHg s / ml -> small enough)
     * double rev;      [1] // extra-splanchnic venous resistance (nominal)
     * double rsv;      [2] // splanchnic venous resistance (nominal)
     * double rmv;      [3] // skeletal muscle venous resistance (nominal)
     * double rhv;      [4] // coronary venous resistance (nominal)
     * double rbv;      [5] // brain venous resistance (nominal)
     * double rpp;      [6] // pulmonary peripheral resistance
     * double rpv;      [7] // pulmonary venous resistance
     * double rsp0;     [8] // extra-splanchnic peripheral resistance    // controlled in cardiovascularControl 
     * double rep0;     [9] // splanchnic peripheral resistance          // controlled in cardiovascularControl
     * double rmp0;     [10] // skeletal muscle peripheral resistance     // controlled in cardiovascularControl 
     * double rhp0;     [11] // coronary peripheral resistance           // controlled in cardiovascularControl 
     * double rbp0;     [12] // brain peripheral resistance               // controlled in cardiovascularControl 
     * double Rpp;      [13] // resistance from pulmonary artery to pulmonary pheripheral compartment (fixed at 0.0001 mmHg s / ml -> small enough)
     * double Rps;      [14] // resistance from pulmonary artery to pulmonary shunt  (fixed at 0.1 mmHg s / ml -> small enough)
     */
    vector<double> Rparam = vector<double>(15);     

    /**
      * @brief Vector of 11 components that contains thoracic veins parameters
      * 
      * double tvd1;    [0]
      * double tvk1;    [1]
      * double tvkr;    [2]
      * double tvvmax;  [3]
      * double tvvu;    [4]
      * double tvr0;    [5]   // nominal thoracic vein resistance
      * double tvd2;    [6]
      * double tvk2;    [7]
      * double tvvmin;  [8]
      * double tvkxp;   [9]
      * double tvkxv;   [10]
    */
    vector<double> TveinParam = vector<double>(11);   

    /**
     * @brief Vector of 15 components that contains heart atria parameters
     * 
     * double vula;     [0]   // left atrium unstressed volume
     * double vura;     [1]   // right atrium unstressed volume
     * double s;        [2]   
     * double trRA;     [3]
     * double trLA;     [4]
     * double tcRA;     [5]
     * double tcLA;     [6]
     * double TrpRA;    [7]
     * double TrpLA;    [8]
     * double TcpRA;    [9]
     * double TcpLA;    [10]
     * double eALA;     [11]
     * double eBLA;     [12]
     * double eARA;     [13]
     * double eBRA;     [14]
    */
    vector<double> atriaParam = vector<double>(15);  

    /**
     * @brief Vector of 10 components that contains heart ventricles parameters
     * 
     * double p0lv;     [0]
     * double p0rv;     [1]
     * double kelv;     [2]
     * double kerv;     [3]
     * double vulv;     [4]   // left ventricle unstressed volume
     * double vurv;     [5]   // right ventricle unstressed volume
     * double emaxlv0;  [6]   // controlled in cardiovascularControl 
     * double emaxrv0;  [7]   // controlled in cardiovascularControl
     * double krlv;     [8]
     * double krrv;     [9]
    */
    vector<double> ventrParam = vector<double>(10);   

    /**
     * @brief Vector of 34 components that contains heart valve parameters
     * 
     * double AmaxTV;        [0]    // maximum area of tricuspid valve
     * double AminTV;        [1]    // minimum area of tricuspid valve
     * double lTV;           [2]    // lenght of tricuspid valve
     * double deltap_openTV; [3]
     * double k_openTV;      [4]
     * double deltap_closeTV;[5]
     * double k_closeTV;     [6]
     *
     * double AmaxPV;        [7]    // maximum area of pulmonary valve
     * double AminPV;        [8]    // minimum area of pulmonary valve
     * double lPV;           [9]    // lenght of pulmonary valve
     * double deltap_openPV; [10]
     * double k_openPV;      [11]
     * double deltap_closePV;[12]
     * double k_closePV;     [13]
     *
     * double AmaxMV;        [14]   // maximum area of mitral valve
     * double AminMV;        [15]   // minimum area of mitral valve
     * double lMV;           [16]   // lenght of mitral valve
     * double deltap_openMV; [17]
     * double k_openMV;      [18]
     * double deltap_closeMV;[19]
     * double k_closeMV;     [20]
     *
     * double AmaxAV;        [21]   // maximum area of aortic valve
     * double AminAV;        [22]   // minimum area of aortic valve
     * double lAV;           [23]   // lenght of aortic valce
     * double deltap_openAV; [24]
     * double k_openAV;      [25]
     * double deltap_closeAV;[26]
     * double k_closeAV;     [27]
     *
     * double rpprox;        [28]
     * double rsprox;        [29]
     * 
     * double ATV;           [30]   // reference area of tricuspid valve
     * double APV;           [31]   // reference area of pulmonary valve
     * double AMV;           [32]   // reference area of mitral valve
     * double AAV;           [33]   // reference area of aortic valve
     */
    vector<double> valveParam = vector<double>(34);     
    
    /**
      * @brief Vector of 6 components that contains additional parameters
      * 
      * double sh;      [0] // pulmonary shunt (%)
      * double rho;     [1] // blood density (g/mL)
      * 
      * Additional parameters for model without regulation
      * double T;       [2] // heart period (if controls are turned off)
      * double tsys0;   [3] // period of systole LV
      * double ksys;    [4] // s^2 constant for determining hr T
      * double pl;      [5] // pleural pressure
    */
    vector<double> params = vector<double>(6);


    //-----------------------------------------------------------------------------------
    // Algebraic variables
    //-----------------------------------------------------------------------------------

    /**
     * @brief Vector of 33 components that containg the systemic circulations algeraic equations
     * 
     * Systemic arteries
     * double psa;       [0]
     * double qsa;       [1]
     * 
     * Systemic peripheral compartments
     * double pep;       [2]
     * double psp;       [3]
     * double pmp;       [4]
     * double php;       [5]
     * double pbp;       [6]
     * double qep;       [7]
     * double qsp;       [8]
     * double qmp;       [9]
     * double qhp;       [10]
     * double qbp;       [11]
     * 
     * Thoracic vein
     * double ptvtm;     [12]
     * double ptv;       [13]
     * double rtv;       [14]
     * 
     * Systemic veins
     * double pev;       [15]
     * double psv;       [16]
     * double pmv;       [17]
     * double phv;       [18]
     * double pbv;       [19]
     * double qev;       [20]
     * double qsv;       [21]
     * double qmv;       [22]
     * double qhv;       [23]
     * double qbv;       [24]
     * double qtv;       [25]
     * double vsv;       [26]
     * double vmv;       [27]
     * 
     * From systemic artieries to peripheral compartments
     * double qepin;     [28]
     * double qspin;     [29]
     * double qmpin;     [30]
     * double qhpin;     [31]
     * double qbpin;     [32]
    */
    vector<double> sysCircAlg = vector<double>(33);

    /**
     * @brief Vector of 10 components that contains the pulmonary circulation algeraic equations
     * 
     * double ppa;       [0]  // pressure in pulmonary arteries
     * double ppv;       [1]  // pressure in pulmonary veins
     * double qpp;       [2]  // flow in pulmonary pheripheral compartments
     * double qps;       [3]  // flow in pulmonary shunt
     * double qpv;       [4]  // flow in pulmonary veins
     * double qpa;       [5]  // flow in pulmonary arteries
     * double ppp;       [6]  // pressure in pulmonary pheripheral compartments
     * double pps;       [7]  // pressure in pulmonary shunt
     * double qppin;     [8]  // flow entering the pulmonary pheripheral compartments
     * double qpsin;     [9]  // flow entering the pulmonary shunt
     */
    vector<double> pulCircAlg = vector<double>(10);

    /**
     * @brief Vector of 6 components that contains the heart valve albegraic equations
     * 
     * double ppprox;    [0]
     * double psprox;    [1]
     * double deltapTV;  [2]
     * double deltapPV;  [3]
     * double deltapMV;  [4]
     * double deltapAV;  [5]
     */
    vector<double> valveAlg = vector<double>(6);

    /**
     * @brief Vector of 13 equations that contains the heart chambers algebraic equations
     * double tsys;     [0] // duration of systole
     * double u;        [1] // parameter representing the fraction of cardiac cyle
     * double phi;      [2] // ventricle activation function
     * 
     * Ventricles
     * double pmaxlv;   [3] // isometric left ventricle  pressure (w/o viscous losses)
     * double pmaxrv;   [4] // isometric right ventricle  pressure (w/o viscous losses)
     * double plv;      [5] // pressure in left ventricle
     * double prv;      [6] // pressure in right ventricle
     * 
     * Atria
     * double pra;      [7] // pressure in right atria
     * double pla;      [8] // pressure in left atria
     * double eLA;      [9]
     * double eRA;      [10]
     * 
     * Ventricles
     * double rlv;      [11] // left ventricle resistance
     * double rrv;      [12] // right ventricle resistance
     */
    vector<double> heartAlg = vector<double>(13);

    //-----------------------------------------------------------------------------------
    // Controlled Variables
    //-----------------------------------------------------------------------------------
    
    /**
      * @brief Vector of 11 components that contains variables controlled by the cardiovascular control
      * double T;         [0]
      * double vuev;      [1]
      * double vusv;      [2]
      * double vumv;      [3]
      * double rep;       [4]
      * double rsp;       [5]
      * double rmp;       [6]
      * double rhp;       [7]
      * double rbp;       [8]
      * double emaxlv;    [9]
      * double emaxrv;    [10]
    */
    vector<double> controlledVars = vector<double>(11);

    //-----------------------------------------------------------------------------------
    // Outputs
    //-----------------------------------------------------------------------------------

    ofstream sampleCVSstateVars;
    ofstream sampleCVSsysCirc;
    ofstream sampleCVSpulCirc;
    ofstream sampleCVSheart;
    ofstream sampleCVSvalve;
    string nameCVS;

    //-----------------------------------------------------------------------------------
    // Methods
    //-----------------------------------------------------------------------------------

    void init(string ifile, string ifileC, string outDir);
    void readCVSparameters(string _file, string _fileC);
    void printCVSparameters();
    void InitialCond(string _file);
    double getqValvedt(double deltap, double Amax, double Amin, double q, double l, double xi);
    double getxiValvesdt(double deltap,double deltap_open, double k_open, double deltap_close,double k_close,double xi);
    double geteAtria(double t, double tr, double Trp,double tc,double Tcp);
    void getAlgebraicRelations();
    void getTimeDerivative();
    void output();
};
