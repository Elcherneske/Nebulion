// Copyright 2023 Jimmy Eng
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#ifndef _COMMON_H_
#define _COMMON_H_


#include "ThreadPool.h"
#include "MSReader.h"
#include "Spectrum.h"
#include "MSObject.h"

/*#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <vector>
#include <utility>
#include <set>
#include <cfloat>
#include <iostream>
#include <map>
#include <algorithm>
#include <numeric>
#include <unordered_map>
#include <random>*/ 


#define nebulion_version   "2024.12 rev. 0"
#define copyright "(c) Peking University"

// Redefined how the bin offset is interpreted and applied.  The valid range for the offset is
// now between 0.0 and 1.0 and scales to the binWidth.
//#define BIN(dMass) (int)((dMass)*Param::g_staticParams.dInverseBinWidth + Param::g_staticParams.dOneMinusBinOffset)

#define isEqual(x, y) (std::abs(x-y) <= ( (std::abs(x) > std::abs(y) ? std::abs(y) : std::abs(x)) * FLT_EPSILON))

#define cometbitset(byte, nbit)   ((byte) |=  (1<<(nbit)))  // https://www.codementor.io/@hbendali/c-c-macro-bit-operations-ztrat0et6
#define cometbitclear(byte, nbit) ((byte) &= ~(1<<(nbit)))
#define cometbitflip(byte, nbit)  ((byte) ^=  (1<<(nbit)))
#define cometbitcheck(byte, nbit) ((byte) &   (1<<(nbit)))

// A common used class for static const value

struct ElementMassUtils
{
   static constexpr double Hydrogen_Mono = 1.007825032;
   static constexpr double Carbon_Mono = 12.0000000000;
   static constexpr double Nitrogen_Mono = 14.0030740044;
   static constexpr double Oxygen_Mono = 15.9949146196;
   static constexpr double Sulfur_Mono = 31.9720711744;
   static constexpr double Phosphorus_Mono = 30.9737619984;
   static constexpr double Selenium_Mono = 79.9165218;
   static constexpr double Sodium_Mono = 22.98976928;
   static constexpr double Potassium_Mono = 39.0983;
   static constexpr double Hydrogen_2_Mono = 2.0141017781;
   static constexpr double Carbon_13_Mono = 13.0033548351;
   static constexpr double Nitrogen_15_Mono = 15.0001088989;
   static constexpr double Oxygen_15_Mono = 17.9991596128;
   static constexpr double Hydrogen_Ave = 1.00794;
   static constexpr double Carbon_Ave = 12.0107;
   static constexpr double Nitrogen_Ave = 14.0067;
   static constexpr double Oxygen_Ave = 15.9994;
   static constexpr double Sulfur_Ave = 32.065;
   static constexpr double Phosphorus_Ave = 30.973761;
   static constexpr double Selenium_Ave = 78.96;

   static constexpr double Proton_Mono = 1.00727646688; //without electron
};

struct CommonValues
{
   static constexpr double dFloatZero = 1e-6;                        // 0.000001
   static constexpr int    iMaxNameSize = 128;                       // max # for every place for name

   //Param for Nucleic Acid
   static constexpr int iMaxNucleicAcidLenght = 512;                 // max # of Nucleotide for a Nucleic Acid

   //Param for Oligonucleotide
   static constexpr int iMinOligonucleotideLength = 1;               // min # of Nucleotide for a Oligonucleotide
   static constexpr int iMaxOligonucleotideLength = 51;              // max # of Nucleotide for a Oligonucleotide
   static constexpr int iAsciiValueSize = 128;                       // store the mass of NT
   static constexpr int iMaxOligonucleotideCharge = 10;              // max charge of Oligonucleotide
   static constexpr double dMaxOligonucleotideMass = 8000.0;         // max mass of Oligonucleotide
   static constexpr double dMinOligonucleotideMass = 600.0;          // min mass of Oligonucleotide

   //Param for Enzyme

   //Param for Fragments
   static constexpr int iMaxFragmentCharge = 4;                      // max charge of fragment
   static constexpr double dMaxFragmentMass = 4000.0;                // max mass of Fragment
   static constexpr double dMinFragmentMass = 200.0;                 // min mass of Fragment

   //File 
   static constexpr int iFileBufferSize = 8192;                      // use for file stream output
   static constexpr int iStreamBufferSize = 2048;                    // for printout stream or other temp char array
   static constexpr int iMaxFileNameSize = 4096;                     // for char array storing filename


   //thread
   static constexpr int iMaxThreadProc = 128;
};

using namespace MSToolkit;
using namespace std;

#endif // _COMMON_H_
