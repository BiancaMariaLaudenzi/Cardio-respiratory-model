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
#include <iostream>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
//#include <experimental/filesystem>
#include <cstdlib>
//#include <direct.h> 

#include "closedloop.h"


int main(int argc, char *argv[])
{ 
        closedloop cl; // declare closedloop of type cl

        // assign values to params not read from file
        int verbose = 1;
        string ifile = argv[1];
                
        // initialise systems (CVS, lungs, gas transport)
        cl.init(ifile, verbose);
        // update systems
        cl.timeLoop(verbose);

        cout << "#######################" << endl;

    return 0;
}
