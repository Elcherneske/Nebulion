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

#include "ModificationUtils.h"
#include "CalculateUtils.h"
#include "Threading.h"
#include "WriteOutUtils.h"

unordered_map<string, int>  ModificationUtils::varModificationMap;
Mutex                       ModificationUtils::_mapAccessMutex;

//static function
static void generateCartesianCombinations(const vector<vector<int>>& lists, vector<int>& currentCombination, vector<vector<int>>& result, int depth) 
{
    if (depth == lists.size())
    {
        result.push_back(currentCombination);
        return;
    }

    for (int i = 0; i < lists[depth].size(); ++i)
    {
        currentCombination[depth] = lists[depth][i];
        generateCartesianCombinations(lists, currentCombination, result, depth + 1);
    }
}

static vector<vector<int>> getAllCombinationSets(int n, int k)
{
    vector<vector<int>> allSets;

    if (n <= 0 || k <= 0)
        return allSets;

    if (k > n)
        k = n;

    for (int i = k; i > 0; i--)
    {
        int nck = CombinatoricsUtils::nChooseK(n, i);

        int** combinationSets = CombinatoricsUtils::makeCombinations(n, i, nck);
        for (int j = 0; j < nck; ++j)
        {
            int* combination = combinationSets[j];
            vector<int> set;
            for (int k = 0; k < i; ++k)
            {
                set.push_back(combination[k]);
            }
            allSets.push_back(set);
        }

        for (int j = 0; j < nck; ++j)
        {
            delete[] combinationSets[j];
        }
        delete[] combinationSets;
    }

    return allSets;
}

//ModificationUtils Function
bool ModificationUtils::AddStaticMod(vector<NucleicAcid>& vNucleicAcidList, string sStaticModFilename)
{
    ifstream inputFile(sStaticModFilename);
    string sLine;

    if(!inputFile)
    {
        string sErrorMsg = "Error: not valid static modification filename: " + sStaticModFilename;
        logerr(sErrorMsg.c_str());
        return false;
    }
    else
    {
        cout << "Begin add static modifications" << endl;
    }

    unordered_map<string, NucleicAcid*> dictionary;

    while (getline(inputFile, sLine)) 
    {
        if (sLine.find('#') != string::npos)
        {
            sLine.erase(sLine.find('#'));
            sLine.erase(std::remove_if(sLine.begin(), sLine.end(), [](char c) {return std::isspace(static_cast<char>(c));}), sLine.end()); //去除两端多余的空格
            if (sLine.size() == 0)
                continue;
        }
        stringstream lineStream(sLine);
        string sNucleicAcidName, sModName, sPosition, sMass;

        lineStream >> sNucleicAcidName >> sModName >> sPosition >> sMass;

        if (dictionary.find(sNucleicAcidName) != dictionary.end())
        {
            if ((*dictionary[sNucleicAcidName]).sName == sNucleicAcidName)
            {
                Modification modification = Modification();
                modification.dModMass = stod(sMass);
                modification.iModID = 0;                    //represent static modification
                modification.sModName = sModName;
                (*dictionary[sNucleicAcidName]).vMods.push_back({stoi(sPosition), modification});
            }
            else
            {
                string sErrorMsg = "Error: find nuclic acid in dictionary error for: " + sNucleicAcidName + "\n";
                logout(sErrorMsg.c_str());
                return false;
            }
        }
        else
        {
            bool bFindTmp = false;
            for (NucleicAcid& nucleicAcid : vNucleicAcidList)
            {
                if (nucleicAcid.sName == sNucleicAcidName)
                {
                    bFindTmp = true;
                    Modification modification = Modification();
                    modification.dModMass = stod(sMass);
                    modification.iModID = 0;                    //represent static modification
                    modification.sModName = sModName;
                    nucleicAcid.vMods.push_back({stoi(sPosition), modification});
                    dictionary[sNucleicAcidName] = &nucleicAcid;
                }
            }

            if (bFindTmp == false)
            {
                string sErrorMsg = "Error: cannot find nuclic acid: " + sNucleicAcidName + " for adding static modification " + sModName + "\n";
                logerr(sErrorMsg.c_str());
                return false;
            }
        }
    }

    cout << "finish adding static modification" << endl;
    return true;
}

void ModificationUtils::AddVarModToOligo(Oligonucleotide& oligonucleotide, vector<VarModInputItem> vModificationList) // generate varmod for single oligo, 这里oligo不需要加锁，因为不同线程用不同oligo
{
    pair<string, int> pair = ModificationUtils::ExtVarModNT(oligonucleotide.sSequence, vModificationList);
    string varSubSequnce = pair.first;
    int iMaxVarModification = pair.second;

    if (iMaxVarModification < 0)
        iMaxVarModification = 0;
    if (iMaxVarModification > Param::g_staticParams.options.iMaxVarModPerOligonucleoitde)
        iMaxVarModification = Param::g_staticParams.options.iMaxVarModPerOligonucleoitde;

    Threading::LockMutex(ModificationUtils::_mapAccessMutex);
    if (!(ModificationUtils::varModificationMap.find(varSubSequnce) == ModificationUtils::varModificationMap.end()))
    {
        int index = (ModificationUtils::varModificationMap.find(varSubSequnce))->second;
        oligonucleotide.iVarModificationIndex = index;
        Threading::UnlockMutex(ModificationUtils::_mapAccessMutex);
    }
    else
    {
        Threading::UnlockMutex(ModificationUtils::_mapAccessMutex);
        // form various modification structure
        VarModification varmodification = VarModification();
        varmodification.sSequence = varSubSequnce;
        varmodification.iMaxModificationNumber = iMaxVarModification;
        ModificationUtils::FormVarMod(varmodification, vModificationList);   //生成varmod的过程可以多线程计算

        //add to Param::g_vVarModification and update map，需要mutex
        Threading::LockMutex(Param::g_vVarModificationMutex);
        Threading::LockMutex(ModificationUtils::_mapAccessMutex);
        if (ModificationUtils::varModificationMap.find(varSubSequnce) == ModificationUtils::varModificationMap.end())
        {
            int index = Param::g_vVarModification.size();
            Param::g_vVarModification.push_back(varmodification);
            ModificationUtils::varModificationMap[varSubSequnce] = index;
            oligonucleotide.iVarModificationIndex = index;
        }
        else
        {
            int index = (ModificationUtils::varModificationMap.find(varSubSequnce))->second;
            oligonucleotide.iVarModificationIndex = index;
        }
        Threading::UnlockMutex(ModificationUtils::_mapAccessMutex);
        Threading::UnlockMutex(Param::g_vVarModificationMutex);
    }
    
    return;
}

bool ModificationUtils::AddVarMod(ThreadPool* _tp, vector<Oligonucleotide>& vOligonucleotideList, string sVarModFilename) // read VarModFile and generate Mod
{
    ThreadPool* pSearchThreadPool = _tp;

    vector<VarModInputItem> vModificationList;
    if(!ModificationUtils::ReadVarModInput(sVarModFilename, vModificationList))
    {
        string sErrorMsg = "Error: read various modification file fail";
        logerr(sErrorMsg.c_str());
        return false;
    }
    
    //先创建varmodification，push到g_varmodification中,这个过程目前不用多线程
    for (int i = 0; i < vOligonucleotideList.size(); i++)
    {
        Oligonucleotide& oligo = vOligonucleotideList[i];
        pair<string, int> pair = ModificationUtils::ExtVarModNT(oligo.sSequence, vModificationList);
        string varSubSequnce = pair.first;
        int iMaxVarModification = pair.second;
        if (iMaxVarModification < 0)
            iMaxVarModification = 0;
        if (iMaxVarModification > Param::g_staticParams.options.iMaxVarModPerOligonucleoitde)
            iMaxVarModification = Param::g_staticParams.options.iMaxVarModPerOligonucleoitde;

        if (!(ModificationUtils::varModificationMap.find(varSubSequnce) == ModificationUtils::varModificationMap.end()))//需要测试一下这个find的速度
        {
            int index = (ModificationUtils::varModificationMap.find(varSubSequnce))->second;
            oligo.iVarModificationIndex = index;
        }
        else
        {
            // form various modification structure
            VarModification varmodification = VarModification();
            varmodification.sSequence = varSubSequnce;
            varmodification.iMaxModificationNumber = iMaxVarModification;
            //这里先不进行 FormVarMod
            //ModificationUtils::FormVarMod(varmodification, vModificationList);   //生成varmod的过程可以多线程计算
            int index = Param::g_vVarModification.size();
            Param::g_vVarModification.push_back(varmodification);
            ModificationUtils::varModificationMap[varSubSequnce] = index;
            oligo.iVarModificationIndex = index;
        }

    }

    //多线程form Varmodification的内容
    for (int i = 0; i < Param::g_vVarModification.size(); i++)
    {
        VarModification& varMod = Param::g_vVarModification[i];
        while (pSearchThreadPool->jobs_.size() >= CommonValues::iMaxThreadProc)
        {
            pSearchThreadPool->wait_for_available_thread();
        }

        pSearchThreadPool->doJob(std::bind(ModificationUtils::FormVarMod, std::ref(varMod), vModificationList));

        //ModificationUtils::AddVarModToOligo(oligo, vModificationList);
        if (!(i % 1000))
        {
            WriteOutUtils::ShowProgressBar(i, Param::g_vVarModification.size(), "Form Variable Modification", "");
        }
    }

    pSearchThreadPool->wait_on_threads();
    WriteOutUtils::ShowProgressBar(Param::g_vVarModification.size(), Param::g_vVarModification.size(), "Form Variable Modification", "done\n");

    Threading::DestroyMutex(ModificationUtils::_mapAccessMutex);
    return true;

}

bool ModificationUtils::FormVarMod(VarModification& varModification, vector<VarModInputItem> vModificationList)
{   
    string sSequence = varModification.sSequence;                       //even sSequence == "", it should also work
    int iMaxModificationNumber = varModification.iMaxModificationNumber;
    if(varModification.iMaxModificationNumber > varModification.sSequence.size())
        varModification.iMaxModificationNumber = varModification.sSequence.size();

    //add one item to represent no various modification,确保至少有个无修饰的项
    vector<pair<int, Modification>> emptyModification;  
    varModification.varModificationList.push_back(emptyModification);

    vector<vector<int>> vvAllCombinationSet = getAllCombinationSets(sSequence.size(), iMaxModificationNumber);
    for (auto& vModPos : vvAllCombinationSet)
    {
        if(!ModificationUtils::FormSingleCombination(vModPos, sSequence, vModificationList, varModification.varModificationList))
        {
            return false;
        }
    }

    for(auto& varModCombination: varModification.varModificationList)
    {
        double dMassSum = 0.0;
        for (auto& pair: varModCombination)
        {
            dMassSum += pair.second.dModMass;
        }
        varModification.vVarModificationMassList.push_back(dMassSum);
    }
    varModification.CalcMassRange();
}

bool ModificationUtils::FormSingleCombination(vector<int> vModPos, string sSequence, vector<VarModInputItem> vModificationList, vector<vector<pair<int, Modification>>>& result)
{
    if (vModPos.size() == 0 || vModificationList.size() == 0)
        return true;

    int iPosLen = vModPos.size();
    
    vector<vector<Modification>> vvModificationSet;     // pos_len * candidate_size
    vector<vector<int>> vvModificationIndexSet;         // pos_len * candidate_size
    vector<vector<int>> vvSelectedCombinations;         // combination_size * pos_len

    for (int& pos: vModPos)
    {
        vector<Modification> vModCandidates;
        vector<int> vCandidateIndexSet;
        char c = sSequence[pos];
        for (auto& varModEntry : vModificationList)
        {
            if (varModEntry.sAllowNT.find(c) != string::npos)
            {
                vModCandidates.push_back(varModEntry.modification);
                vCandidateIndexSet.push_back(vCandidateIndexSet.size());
            }
        }
        vvModificationIndexSet.push_back(vCandidateIndexSet);
        vvModificationSet.push_back(vModCandidates);
    }

    vector<int> currentCombination(iPosLen);
    generateCartesianCombinations(vvModificationIndexSet, currentCombination, vvSelectedCombinations, 0);

    for (auto& vSelectedCombination : vvSelectedCombinations)
    {
        vector<pair<int, Modification>> vModificationCombination;
        for (int i = 0; i < vSelectedCombination.size(); i++)
        {
            int modificationIndex = vSelectedCombination[i];
            int NTPosition = vModPos[i];
            vModificationCombination.push_back({NTPosition, vvModificationSet[i][modificationIndex]});
        }
        result.push_back(vModificationCombination);
    }

    return true;
}

bool ModificationUtils::ReadVarModInput(string sFilename, vector<VarModInputItem>& vModificationList)
{
    ifstream inputFile(sFilename);
    string sLine;
    int iModID = 0;

    if(!inputFile)
    {
        string sErrorMsg = "Error: not valid various modification filename: " + sFilename;
        logerr(sErrorMsg.c_str());
        return false;
    }

    unordered_map<string, int> duplicateDictionary;             //used to reduce duplicate

    while (getline(inputFile, sLine)) 
    {
        if (sLine.find('#') != string::npos)
        {
            sLine.erase(sLine.find('#'));
            sLine.erase(std::remove_if(sLine.begin(), sLine.end(), [](char c) {return std::isspace(static_cast<char>(c));}), sLine.end()); //去除两端多余的空格
            if (sLine.size() == 0)
                continue;
        }
        stringstream lineStream(sLine);
        string sModName, sMass, sMaxNum, sAllowNT;

        lineStream >> sModName >> sMass >> sAllowNT >> sMaxNum;

        if (duplicateDictionary.find(sModName) != duplicateDictionary.end())
        {
            string sErrorMsg = "Error: duplicate various modification: " + sModName + "\n" + "ignore this duplicate various modification.\n";
            logerr(sErrorMsg.c_str());
        }
        else
        {
            duplicateDictionary[sModName] = 1;                  // mark this modification

            Modification modification = Modification();
            modification.sModName = sModName;
            modification.dModMass = stod(sMass);
            modification.iModID = ++iModID;

            int iMaxNum = stoi(sMaxNum);
            if (iMaxNum < 0)
                iMaxNum = 0;
            if (iMaxNum > Param::g_staticParams.options.iMaxVarModPerOligonucleoitde)
                iMaxNum = Param::g_staticParams.options.iMaxVarModPerOligonucleoitde;

            VarModInputItem item;
            item.sAllowNT = sAllowNT;
            item.modification = modification;
            item.iMaxNum = iMaxNum;

            vModificationList.push_back(item);
        }
    }

    return true;
}

pair<string, int> ModificationUtils::ExtVarModNT(string sSequence, vector<VarModInputItem> vModificationList)
{
    string varSequence = "";
    int iMaxNum = 0;

    for (int i = 0; i < sSequence.size(); i++)
    {
        char NT = sSequence[i];
        bool bPush = false;
        for (auto& mod : vModificationList)
        {
            if (mod.sAllowNT.find(NT) != string::npos)
            {
                bPush = true;
                if (mod.iMaxNum > iMaxNum)
                    iMaxNum = mod.iMaxNum;
            }
        }

        if (bPush)
            varSequence.push_back(NT);
    }

    return {varSequence, iMaxNum};
}
