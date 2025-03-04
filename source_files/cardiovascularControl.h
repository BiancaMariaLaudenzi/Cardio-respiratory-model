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
 * @brief Cardiovascular control model
 * 
 * @note Cardiovascular control system models proposed in: 
 * 
 * - M. Ursino and E. Magosso. Acute cardiovascular response to isocapnic hypoxia. I. A mathematical
 *   model. American Journal of Physiology. Heart and Circulatory Physiology, 279(1):H149–165, July 2000.
 * 
 * - E. Magosso and M. Ursino. A mathematical model of CO2 effect on cardiovascular regulation.
 *   American Journal of Physiology. Heart and Circulatory Physiology, 281(5):H2036–2052, November 2001.
 * 
 * - Mauro Ursino and Elisa Magosso. A theoretical analysis of the carotid body chemoreceptor
 *    response to O2 and CO2 pressure changes. Respiratory Physiology & Neurobiology, 130(1):99–110, March 2002.
 */
class cardiovascularControl
{
    public:
        /// Time in the lungs
        double time;
        /// dt used for the update
        double dt_control;
        /// dt used for sampling
        double dtSample;
        /// time at which sampling starts
        double tSampleIni;
        /// time at which sampling ends
        double tSampleEnd;
        /// time at which simulation starts
        double tIni;
        /// iteration number
        double iT;
        /// Integer determining whether outputs are printed to teminal/log file
        int verbose;

        //-----------------------------------------------------------------------------------
        // Parameters in common wih other classes
        //-----------------------------------------------------------------------------------

        /// @brief  basal heart period
        double T0; 
        /**
         * @brief Common params 
         * rmp0 [0] 
         * rsp0 [1]
         * rep0 [2]
         * vumv0 [3]
         * vusv0 [4]
         * vuev0 [5]
         * emaxlv0 [6]
         * emaxrv0 [7]
         * rbp0 [8]
         * rhp0 [9]
         */
        vector <double> commonParams = vector <double> (10);

        /**
         * @brief Common control params
         * - double paco2n [0] 
         * - double caco2n [1]
         */
        vector <double> commonControlParams = vector <double> (2);

        //-----------------------------------------------------------------------------------
        // Local metabolic regulation
        //-----------------------------------------------------------------------------------

        /**
         * @brief Brain autoregulation parameters
         * - double gbpn [0]
         * - double gbo2 [1]
         * - double cvbo2n [2] 
         * - double A [3]
         * - double B [4]
         * - double CC [5]
         * - double D [6]
         */
        vector <double> brainAutoParams = vector <double> (7);

        /**
         * @brief  Coronary autoregulation params
         * - double gho2 [0] 
         * - double gmo2 [1] 
         * - double cvho2n [2]
         * - double cvmo2n [3] 
         * - double khco2 [4]
         * - double kmco2 [5]
         */
        vector <double> coromuscleAutoParams = vector <double> (6);

        /**
         * @brief Common autoregulation params
         * - double tauo2A [0]
         * - double tauco2A [1]
         */
        vector <double> commonAutoParams = vector <double> (2);

        //-----------------------------------------------------------------------------------
        // Afferent pathways
        //-----------------------------------------------------------------------------------

        /**
         * @brief Vector that contains parameters related to the afferent baroreflex pathway
         * - double tauzB [0] 
         * - double taupB [1] 
         * - double fabmin [2] 
         * - double fabmax [3] 
         * - double pn [4] 
         * - double kab [5]
         */
        vector <double> baroAfferentParams = vector <double> (6);
        
        /**
         * @brief Vector that contains parameters related to the chemoreflex pathway 
         * - double Ac [0]
         * - double Bc [1]
         * - double kco2 [2]
         * - double ko2 [3]
         * - double kstat [4]
         * - double kdyn [5]
         * - double tauCAP [6]
         * - double tauccCO2dyn [7] 
         */
        vector <double> chemoAfferentParams = vector <double> (8);

        /**
         * @brief vector containing lung stretch receptor parameters
         * - double gap [0]
         * - double tauLSR [1]
         */
        vector <double> lsrAfferentParams = vector <double> (2);

        //-----------------------------------------------------------------------------------
        // Efferent pathways
        //-----------------------------------------------------------------------------------
        // h = heart, p = peripheral, v = veins

        /**
         * @brief Parameters for firing rates of efferent sympathetic pathway
         * - double fesinf [0]
         * - double fes0 [1]
         * - double fesmax [2]
         * - double kes [3]
         * - double Wbsh [4]
         * - double Wbsp [5]
         * - double Wbsv [6]
         * - double Wcsh [7]
         * - double Wcsp [8]
         * - double Wcsv [9]
         * - double Wpsh [10]
         * - double Wpsp [11]
         * - double Wpsv [12]
         */
        vector <double> efferentSympParams = vector <double> (13);

        /**
         * @brief Parameters for firing rates of efferent vagal pathway
         * - double fevinf [0]
         * - double fev0 [1]
         * - double fab0 [2]
         * - double kev [3]
         * - double Wcv [4]
         * - double Wpv [5]
         * - double thetav [6]
         */
        vector <double> efferentVagalParams = vector <double> (7);


        //-----------------------------------------------------------------------------------
        // CNS ischemic response
        //-----------------------------------------------------------------------------------
        // h = heart, p = peripheral, v = veins

        /**
         * @brief Parameters for CNS ischemic response
         * - double chish [0]
         * - double chisp [1]
         * - double chisv [2]
         * - double Ptildeo2sh [3]
         * - double Ptildeo2sp [4]
         * - double Ptildeo2sv [5]
         * - double kiscsh [6]
         * - double kiscsp [7]
         * - double kiscsv [8]
         * - double tauisc [9]
         * - double taucc [10]
         * - double thetashn [11]
         * - double thetaspn [12]
         * - double thetasvn [13]
         * - double gccsh [14]
         * - double gccsp [15]
         * - double gccsv [16]
         */
        vector <double> CNSResponseParams = vector <double> (17);

        //-----------------------------------------------------------------------------------
        // Reflex effectors
        //-----------------------------------------------------------------------------------
        
        /**
         * @brief Parameters describing changes in rmp, rsp, rep
         * - double Grmp [0]
         * - double Grsp [1]
         * - double Grep [2]
         * - double Drmp [3]
         * - double Drsp [4]
         * - double Drep [5]
         * - double fesmin [6]
         * - double taurmp [7]
         * - double taursp [8]
         * - double taurep [9]
         */
        vector <double> resControlParams = vector <double> (10);

        /**
         * @brief Parameters describing changes in vumv, vusv,vuev
         * - double Gvumv [0]
         * - double Gvusv [1]
         * - double Gvuev [2]
         * - double Dvumv [3]
         * - double Dvusv [4]
         * - double Dvuev [5]
         * - double tauvumv [6]
         * - double tauvusv [7]
         * - double tauvuev [8]
         * note: fesmin was defined for resistances
         */
        vector <double> vuControlParams = vector <double> (9);

        /**
         * @brief Parameters describing changes in cardiac properties
         * - double Gemaxlv [0]
         * - double Gemaxrv [1]
         * - double Demaxlv [2]
         * - double Demaxrv [3]
         * - double tauemaxlv [4]
         * - double tauemaxrv [5]
         * - double GTs [6]
         * - double DTs [7]
         * - double tauTs [8]
         * - double GTv [9]
         * - double DTv [10]
         * - double tauTv [11]
         * note: fesmin was defined for resistances
         */
        vector <double> heartControlParams = vector <double> (12);

        //-----------------------------------------------------------------------------------
        // State Variables
        //-----------------------------------------------------------------------------------

        /**
         * @brief Differential variables: Autoregulation
         * - double xbO2 [0]
         * - double xbCO2 [1]
         * - double xhO2 [2]
         * - double xhCO2 [3]
         * - double xmO2 [4]
         * - double xmCO2 [5]
         */
        vector <double> stateVarAuto = vector <double> (6);
        vector <double> dstateVarAuto = vector <double> (6);

        /**
         * @brief Differential variables: Baroreflex
         * - double Ptilde
         */
        double stateVarBaro;
        double dstateVarBaro;

        /**
         * @brief Differential variables: Chemoreflex 
         * - double phico2Dyn [0]
         * - double phiapc [1]
         */
        vector <double> stateVarChemo = vector <double> (2);
        vector <double> dstateVarChemo = vector <double> (2);

        /**
         * @brief Differential variables: LSR
         * - double fap 
         */
        double stateVarLSR;
        double dstateVarLSR;

        /**
         * @brief Differential variables: CNS ischemic response
         * - double deltathetao2sh [0]
         * - double deltathetaco2sh [1]
         * - double deltathetao2sp [2]
         * - double deltathetaco2sp [3]
         * - double deltathetao2sv [4]
         * - double deltathetaco2sv [5]
         */
        vector <double> stateVarCNS = vector <double> (6);
        vector <double> dstateVarCNS = vector <double> (6);

        /**
         * @brief Differential variables: Sympathetic and parasympathetic changes in resistances, unstressed volumes and cardiac proprerties
         * - double deltarmp [0]
         * - double deltarsp [1]
         * - double deltarep [2]
         * - double deltavumv [3]
         * - double deltavusv [4]
         * - double deltavuev [5]
         * - double deltaemaxlv [6]
         * - double deltaemaxrv [7]
         * - double deltaTs [8]
         * - double deltaTv [9]
         */
        vector <double> stateVarSVResponse = vector <double> (10);
        vector <double> dstateVarSVResponse = vector <double> (10);

        vector <double> stateVars;
        vector <double> dvdt;

        //-----------------------------------------------------------------------------------
        // Algebraic variables
        //-----------------------------------------------------------------------------------

        /**
         * @brief Algebraic relations: Autoregulation
         * - double gbp [0]
         * - double phib [1]
         * - double phih [2]
         * - double phim [3]
         * - double rbp [4]
         * - double rhp [5]
         * - double rmp [6]
         */
        vector <double> algVarAuto = vector <double> (7);

        /**
         * @brief Algebraic relations: Baroreflex
         * - double fab
         */
        double algVarBaro;

        /**
         * @brief Algebraic relations: Chemoreflex
         * - double Xo2 [0]
         * - double phistat [1]
         * - double fcstat [2]
         * - double fcdyn [3]
         * - double fapc [4]
         *          */
        vector <double> algVarChemo = vector <double> (5);

        /**
         * @brief Algebraic relations: LSR
         * - double phiap
         */
        double algVarLSR;

        /**
         * @brief Algebraic relations: firing rates of efferent pathways
         * - double fsh [0]
         * - double fsp [1]
         * - double fsv [2]
         * - double fv [3]
         */
        vector <double> algVarFiring = vector <double> (4);

        /**
         * @brief Algebraic relations: CNS ischemic response
         * - double omegash [0]
         * - double omegasp [1]
         * - double omegasv [2]
         * - double thetash [3]
         * - double thetasp [4]
         * - double thetasv [5]
         */
        vector <double> algVarCNS = vector <double> (6);

        /**
         * @brief Algebraic variables: Sympathetic and parasympathetic changes in resistances, unstressed volumes and cardiac proprerties
         * - double sigmarmp [0]
         * - double sigmarsp [1]
         * - double sigmarep [2]
         * - double sigmavumv [3]
         * - double sigmavusv [4]
         * - double sigmavuev [5]
         * - double sigmaemaxlv [6]
         * - double sigmaemaxrv [7]
         * - double sigmaTs [8]
         * - double sigmaTv [9]
         */ 
        vector <double> algVarSVResponse = vector <double> (10);

<<<<<<< HEAD
        vector <double> algVars = vector <double> (45);
=======
        vector <double> algVars;
>>>>>>> origin/main
        
        //-----------------------------------------------------------------------------------
        // Variables for and from other classes
        //-----------------------------------------------------------------------------------

        /**
         * @brief Variables from transport class
         * - double paco2 [0]
         * - double pao2 [1]
         * - double sao2 [2]
         * - double caco2 [3]
         * - double cbpo2 [4]
         * - double chpo2 [5]
         * - double cmpo2 [6]
         * - double dcaco2dt [7]
         */
        vector <double> varsTransport = vector <double> (8);

        /**
         * @brief Variables from lungMechanics class
         * - double Vt
         */
        double Vt;

        /**
         * @brief Variables from cardiovascular class
         * - double Pb
         * - double dPbdt
         */
        double Pb;
        double dPbdt;

        /**
         * @brief Variables delayed
         * - double fsprmpDel [0]
         * - double fsprspDel [1]
         * - double fsprepDel [2]
         * - double fsvvumvDel [3]
         * - double fsvvusvDel [4]
         * - double fsvvuevDel [5]
         * - double fshemaxlvDel [6]
         * - double fshemaxrvDel [7]
         * - double fshTsDl [8]
         * - fouble fvTvDel [9]
        */
        vector <double> delayedVars = vector <double> (10);

        /**
         * @brief Controlled variables
         * - double rmp [0] 
         * - double rsp [1]
         * - double rep [2]
         * - double vumv [3]
         * - double vusv [4]
         * - double vuev [5]
         * - double emaxlv [6]
         * - double emaxrv [7]
         * - double T [8] 
         * - double rbp [9]
         * - double rhp [10]
         */
        vector <double> controlledVars = vector <double> (11);

        //-----------------------------------------------------------------------------------
        // Methods
        //-----------------------------------------------------------------------------------

        void init(string ifile, string ifileC, string outDir);
<<<<<<< HEAD
        void readCardioControlParameters(string _file, string _fileC);
        void printCardioControlParameters();
        void InitialCond(string _file);
        void getAlgebraicRelationsAuto();
        void getAlgebraicRelationsBaro();
        void getAlgebraicRelationsChemo();
        void getAlgebraicRelationsLSR();
        void getAlgebraicRelationsCNS();
        void getAlgebraicRelationsFiring();
        void getAlgebraicRelationsSVResponse();
        void getAlgebraicRelationscontrolledVars();
        void getAlgebraicRelations();
        void getTimeDerivativeAuto();
        void getTimeDerivativeBaro();
        void getTimeDerivativeChemo();
        void getTimeDerivativeLSR();
        void getTimeDerivativeCNS();
        void getTimeDerivativeSVResponse();
        void getTimeDerivative();
        void setStateVars();
        void setStateVarsVector();
=======

        void readCardioControlParameters(string _file, string _fileC);

        void printCardioControlParameters();
        
        void InitialCond(string _file);

        void getAlgebraicRelationsAuto();

        void getAlgebraicRelationsBaro();

        void getAlgebraicRelationsChemo();

        void getAlgebraicRelationsLSR();

        void getAlgebraicRelationsCNS();

        void getAlgebraicRelationsFiring();

        void getAlgebraicRelationsSVResponse();

        void getAlgebraicRelationscontrolledVars();
        
        void getAlgebraicRelations();

        void getTimeDerivativeAuto();

        void getTimeDerivativeBaro();

        void getTimeDerivativeChemo();

        void getTimeDerivativeLSR();

        void getTimeDerivativeCNS();

        void getTimeDerivativeSVResponse();

        void getTimeDerivative();

        void setStateVars();

        void setStateVarsVector();

>>>>>>> origin/main
        void output();

        ofstream sampleAlgVars;
        ofstream sampleDiffVars;
        ofstream sampleOtherVars;
        ofstream sampleControlVars;
        string nameControl; 


};      