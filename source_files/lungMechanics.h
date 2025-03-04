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
 * @brief Lung mechanics model 
 * 
 * @note  Lung mechanics model proposed in: 
 * 
 * A Albanese, L Cheng, M Ursino & N. W. Chbat.
 * An integrated mathematical model of the human cardiopulmonary system: model development.
 * American Journal of Physiology-Heart and Circulatory Physiology 2016 310, H899-H921.  
 */

class lungMechanics
{
  public:
    /// times
    double time;

    int verbose;

    //-----------------------------------------------------------------------------------
    // State Variables
    //-----------------------------------------------------------------------------------

    /**
     *  @brief Vector that contains volumes
     * - [0] double vl: volume in the larynx
     * - [1] double vtr: volume in the trachea
     * - [2] double vb: volume in the bronchea
     * - [3] double vA: volume in the alveoli
     * - [4] double vpl: volume in the pleural cavity 
     * - [5] xi
     */
    vector <double> vol = vector <double> (6);

    /**
     * @brief Vector that contains time derivative of volumes
     * - [0] double dvldt: dt volume in the larynx
     * - [1] double dvtrdt: dt volume in the traches
     * - [2] double dvbdt: dt volume in the bronchea
     * - [3] double dvAdt: dt volume in the alveoli
     * - [4] double dvpldt: dt volume in the pleural cavity 
     * - [5] dxidt
     */
    vector <double> dvdt = vector <double> (6);
    

    //-----------------------------------------------------------------------------------
    // Parameters
    //-----------------------------------------------------------------------------------

    /**
     * @brief Vector containing the compliances 
     * - [0] double cl: larynx compliance
     * - [1] double ctr: trachea compliance
     * - [2] double cb: bronchea compliance
     * - [3] double cA: alveoli compliance
     * - [4] double ccw: chest wall compliance
     * - [5] double cabd: abdominal compliance
     */
    vector <double> C = vector <double> (6); 
    
    /**
     * @brief Vector containing the unstressed volumes 
     * - [0] double vul: larynx unstressed volume
     * - [1] double vutr: trachea unstressed volume
     * - [2] double vub: bronchea unstressed volume
     * - [3] double vuA: alveoli unstressed volume
     */
    vector <double> Vu = vector <double> (4);

    /**
     * @brief Vector containing the resistances
     * - [0] double rml: resistance mouth to larynx
     * - [1] double rlt: resistance larynx to trachea
     * - [2] double rtb: resistance trachea to bronchea
     * - [3] double rbA: resistance bronchea to alveoli
     */
    vector <double> R = vector <double> (4);

    /**
     * @brief Respiratory parameters
     * - [0] double rr: respiratory rate
     * - [1] double ieratio: inspiratory-expiratory time ratio 
     * - [2] double frc: functional residual capacity
     * - [3] double pplee: pleural pressure value at end expiration
     * - [4] double pmusmin: minimum end inspiratory pressure
     * - [5] double tauCoeff: time constant coefficient
     * - [6] double pao: airway opening pressure
     * - [7] double pvent: ventilator pressure
     * - [8] double patm: atmospheric pressure
     * - [9] double IAPee: end-expiratory IAP
    */
    vector <double> resp = vector <double> (10);
    
    /**
     * @brief Vector with 6 components that contains pressures
     * - [0] double pl: larynx pressure
     * - [1] double ptr: trachea pressure
     * - [2] double pb: bronchea pressure
     * - [3] double pA: alveolar pressure
     * - [4] double ppl1: (1/ccw * vpl) + pplee + pmus
     * - [5] double pabd: abdominal pressure
     */
    vector <double> P = vector <double> (6);

    
    // set time constants
    /// respiratory cycle duration
    double T;
    /// expiratory fraction              
    double te;
    /// inspiratory fraction                
    double ti;
    /// time constant of the exponential expiratory profile                
    double tau;
    /// local time variable that denotes the "position" within the respiratory cycle and is used when defining pmus               
    double timeloc;        
    /// auxiliary time variable
    double u;   

    //-----------------------------------------------------------------------------------
    // Algebraic Variables
    //-----------------------------------------------------------------------------------

    // Respiratory muscles pressure
    double Pmus;             
    /// total airflow
    double Vdot;   
    /// alveolar airflow           
    double VAdot;
    /// dead space volume             
    double vd;    
    /// total volume            
    double v;                 

    //-----------------------------------------------------------------------------------
    // Outputs
    //-----------------------------------------------------------------------------------

    ofstream sampleLungVol;
    ofstream sampleLungPr;
    ofstream sampleLungFlow;
    string nameLung;

    //-----------------------------------------------------------------------------------
    // Methods
    //-----------------------------------------------------------------------------------

    void init(string ifile, string ifileC, string outDir);
    void readLungMechParameters(string _file, string _fileC);
    void InitialCond();
    void printLungMechParameters();
    void getTimeConstants();
    void getPmus();
    void getAlgebraicRelations();
    void getTimeDerivative();
    void output();
}; 