// Copyright 2024 Yucheng Liao
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


#ifndef _MASSUTILS_H_
#define _MASSUTILS_H_

#include "DataStructure.h"

#include <map>
#include <string>
#include <vector>

class MassUtils
{
public:
    MassUtils();
    ~MassUtils();

    double CalOligoMass(const Oligonucleotide& oligo, const OligonucleotideID& oligoID, bool bAddTermMod, bool bAddStaticMod, bool bAddVarMod);

    double CalFragmentMass(const OligoFragment& fragment, bool bAddTermMod, bool bAddMod);

    std::vector<Modification> GenerateEndModification(std::string sCandidates);

    double CalCutModMass(IonSpecies& ionSpecies);
    double CalCutModMass(IonSpecies& ionSpecies, char cBase);

    double MassToMZ(double dMass, int iCharge, std::string sAdductIon);

    double MZToMass(double dMZ, int iCharge, std::string sAdductIon);

    double GetBaseMass(char cBase);

    double GetNTMass(char cBase);

    int CalBinIndex(double dMZ);

    static MassUtils& getInstance() {
        static MassUtils instance;
        return instance;
    }

private:
    double dBackbone;
    double pdNTMass[CommonValues::iAsciiValueSize];
    double pdBaseMass[CommonValues::iAsciiValueSize];

    std::map<std::string, double> mapAdductMass;
};

#endif
