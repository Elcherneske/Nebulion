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

class ModificationUtils
{
    struct VarModInputItem
    {
        string sAllowNT;
        Modification modification;
        int iMaxNum;

        VarModInputItem()
        {
            sAllowNT = "";
            modification = Modification();
            iMaxNum = 0;
        }

        VarModInputItem(const VarModInputItem& other)
        {
            sAllowNT = other.sAllowNT;
            modification = other.modification;
            iMaxNum = other.iMaxNum;
        }

        VarModInputItem& operator=(const VarModInputItem& other)
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

    static bool AddStaticMod(vector<NucleicAcid>& vNucleicAcidList, string sStaticModFilename);

    static void AddVarModToOligo(Oligonucleotide& oligonucleotide, vector<VarModInputItem> vModificationList);

    static bool AddVarMod(ThreadPool* _tp, vector<Oligonucleotide>& vOligonucleotideList, string sVarModFilename);

private:

    static string ExtVarModNT(string sSequence, vector<string> vAllowModificationNT);
    static pair<string, int> ExtVarModNT(string sSequence, vector<VarModInputItem> vModificationList);

    static bool ReadVarModInput(string sFilename, vector<VarModInputItem>& vModificationList);

    static bool ReadModList(string sFilename);

    static bool FormVarMod(VarModification& varModification, vector<VarModInputItem> vModificationList);

    static bool FormSingleCombination(vector<int> vModPos, string sSequence, vector<VarModInputItem> vModificationList, vector<vector<pair<int, Modification>>>& result);

    static unordered_map<string, int> varModificationMap;
    static Mutex _mapAccessMutex;
}; 

#endif