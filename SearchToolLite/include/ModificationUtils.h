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

#ifndef _MODIFICATIONUTILS_H_
#define _MODIFICATIONUTILS_H_

#include "Common.h"
#include "DataStructure.h"
#include "Param.h"

#include <map>
#include <string>

class ModificationUtils
{
    struct VarModEntry
    {
        std::string sAllowNT;
        Modification modification;
        int iMaxNum;

        VarModEntry()
        {
            sAllowNT = "";
            modification = Modification();
            iMaxNum = 0;
        }

        VarModEntry(const VarModEntry& other)
        {
            sAllowNT = other.sAllowNT;
            modification = other.modification;
            iMaxNum = other.iMaxNum;
        }

        VarModEntry& operator=(const VarModEntry& other)
        {
            if (this != &other)
            {
                sAllowNT = other.sAllowNT;
                modification = other.modification;
                iMaxNum = other.iMaxNum;
            }
            return *this;
        }
    };

public:

    bool AddStaticMod(vector<NucleicAcid>& vNucleicAcidList, const std::string& sStaticModFilename);

    bool AddVarMod(ThreadPool* _tp, vector<Oligonucleotide>& vOligo, vector<VarModFormula>& vVarMods, const std::string& sVarModFilename);

    static ModificationUtils& getInstance()
    {
        static ModificationUtils instance;
        return instance;
    }

private:
    std::map<std::string, std::vector<std::vector<int>>> ExtVarModFormula(const string& sSequence, const vector<VarModEntry>& vVarMods);
    // {formula, [positionIndex1, positionIndex2, ...]}

    bool FormVarMod(VarModFormula& varMod, vector<VarModEntry> vVarModEntrys);

    std::map<std::string, int> mVarModIndex;
    Mutex _mapAccessMutex;
}; 

#endif