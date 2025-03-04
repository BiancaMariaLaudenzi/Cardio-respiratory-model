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
 * @brief Respiratory control model
 * 
 * @note Respiratory control system model proposed in: 
 * 
 * A Albanese, L Cheng, M Ursino & N. W. Chbat.
 * An integrated mathematical model of the human cardiopulmonary system: model development.
 * American Journal of Physiology-Heart and Circulatory Physiology 2016 310, H899-H921.  
 */ 

class respiratoryControl
{
    public:
        /// Time in the lungs
        double time;
        double dt_control;
        double dtSample;
        double tSampleIni;
        double tSampleEnd;
        double tIni;
        double iT;
        int verbose;

        //-----------------------------------------------------------------------------------
        // Parameters
        //-----------------------------------------------------------------------------------

        /**
         * @brief parameters
         * - double Dp [0]
         * - double Gpa [1]
         * - double Gpf [2]
         * - double taupa1 [3]
         * - double taupa2 [4]
         * - double taupf1 [5]
         * - double taupf2 [6]
         * - double fapcn [7]
         * - double Dc [8]
         * - double Gca1 [9]
         * - double Gca2 [10]
         * - double Gcf1 [11]
         * - double Gcf2 [12]
         * - double tauca1 [13]
         * - double tauca2 [14]
         * - double taucf1 [15]
         * - double taucf2 [16]
         * - double paco2n [17]
         * - double Pmusmin0 [18]
         * - double RR0 [19]
         */
        vector <double> Params = vector <double> (20);

        //-----------------------------------------------------------------------------------
        // Delayed variables
        //-----------------------------------------------------------------------------------

        /**
         * @brief Delayed variables from other classes
         * - double paco2Del [0]
         * - double fapcDel [1]
         */
        vector <double> delayedVars = vector <double> (2);

        //-----------------------------------------------------------------------------------
        // Algebraic variables
        //-----------------------------------------------------------------------------------

        /**
         * @brief Algebraic variables
         * - double Pmusmin [0]
         * - double RR [1]
         * - double uc [2]
         * - double up [3]
         */
        vector <double> algVars = vector <double> (4);

        //-----------------------------------------------------------------------------------
        // State Variables
        //-----------------------------------------------------------------------------------

        /**
         * @brief State variables
         * - double DeltaPmusminp [0]
         * - double DeltaPmusminc [1]
         * - double DeltaRRp [2]
         * - double DeltaRRc [3]
         */
        vector <double> stateVars = vector <double> (4);

        vector <double> dvdt = vector <double> (4);

        //-----------------------------------------------------------------------------------
        // Methods
        //-----------------------------------------------------------------------------------

        void init(string ifile, string ifileC, string outDir);
        void readRespControlParameters(string _file, string _fileC);
        void InitialCond(string _file);
        void printRespControlParameters();
        void getAlgebraicRelations();
        void getTimeDerivative();
        void output();

        //-----------------------------------------------------------------------------------
        // Outputs
        //-----------------------------------------------------------------------------------

        ofstream sampleAlgVars;
        ofstream sampleStateVars;
        string nameControl; 

};