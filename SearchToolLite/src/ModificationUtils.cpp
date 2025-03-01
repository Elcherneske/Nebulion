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

#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>

//static function
static void cartesianCombinations(const vector<vector<int>>& lists, vector<int>& currentCombination, vector<vector<int>>& result, int depth) 
{
    if (depth == lists.size())
    {
        result.push_back(currentCombination);
        return;
    }

    for (int i = 0; i < lists[depth].size(); ++i)
    {
        currentCombination[depth] = lists[depth][i];
        cartesianCombinations(lists, currentCombination, result, depth + 1);
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
        int nck = CombinatoricsUtils::getInstance().nChooseK(n, i);

        int** combinationSets = CombinatoricsUtils::getInstance().makeCombinations(n, i, nck);
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
bool ModificationUtils::AddStaticMod(vector<NucleicAcid>& vNucleicAcidList, const std::string& sStaticModFilename)
{
    ifstream inputFile(sStaticModFilename);
    string sLine;

    if(!inputFile.is_open())
    {
        std::cerr << "Error: not valid static modification filename: " + sStaticModFilename << std::endl;
        return false;
    }
    else
    {
        std::cout << "Begin add static modifications" << std::endl;
    }

    TSVParser tsvParser;
    if (!tsvParser.readFile(sStaticModFilename))
    {
        std::cerr << "Error: parse static modification file failed" << std::endl;
        return false;
    }

    std::vector<std::map<std::string, std::string>> vData = tsvParser.getData();
    std::map<std::string, int> mNucleicAcidIndex;

    for (const auto& row : vData)
    {
        std::string sNucleicAcidName = row.at("Nucleic_Acid_Name");
        std::string sModName = row.at("Mod_Name");
        int iPosition = std::stoi(row.at("position"));
        double dMass = std::stod(row.at("mass"));

        if (mNucleicAcidIndex.find(sNucleicAcidName) != mNucleicAcidIndex.end())
        {
            NucleicAcid& nucleicAcid = vNucleicAcidList[mNucleicAcidIndex[sNucleicAcidName]];
            if (nucleicAcid.sName != sNucleicAcidName)
                throw std::runtime_error("Error: nucleic acid name mismatch after indexing");

            Modification modification;
            modification.sModName = sModName;
            modification.dModMass = dMass;
            nucleicAcid.vMods.push_back({iPosition, modification});
        }
        else{
            for (size_t index = 0; index < vNucleicAcidList.size(); ++index)
            {
                if (vNucleicAcidList[index].sName == sNucleicAcidName)
                {
                    mNucleicAcidIndex[sNucleicAcidName] = index;
                    Modification modification;
                    modification.sModName = sModName;
                    modification.dModMass = dMass;
                    vNucleicAcidList[index].vMods.push_back({iPosition, modification});
                }
            }
        }
    }

    std::cout << "finish adding static modification" << std::endl;
    return true;
}

bool ModificationUtils::AddVarMod(ThreadPool* _tp, vector<Oligonucleotide>& vOligo, vector<VarModFormula>& vVarMods, const std::string& sVarModFilename) // read VarModFile and generate Mod
{
    // read var mod file and generate var mod entry
    TSVParser tsvParser;
    if (!tsvParser.readFile(sVarModFilename))
    {
        std::cerr << "Error: parse various modification file failed" << std::endl;
        return false;
    }
    std::vector<std::map<std::string, std::string>> vData = tsvParser.getData();
    vector<VarModEntry> vVarModEntrys;
    size_t iModID = 0;
    for (const auto& row : vData)
    {
        std::string sModName = row.at("Mod_Name");
        std::string sMass = row.at("Mod_Mass");
        std::string sAllowNT = row.at("AllowNT");
        std::string sMaxNum = row.at("max_mods_per_Oligo");

        VarModEntry varModEntry;
        varModEntry.sAllowNT = sAllowNT;
        varModEntry.iMaxNum = std::stoi(sMaxNum);
        varModEntry.modification.sModName = sModName;
        varModEntry.modification.dModMass = std::stod(sMass);
        varModEntry.modification.iModID = iModID++;
        vVarModEntrys.push_back(varModEntry);
    }

    // extract var mod formula for each oligo
    for (int i = 0; i < vOligo.size(); i++)
    {
        Oligonucleotide& oligo = vOligo[i];
        std::map<std::string, std::vector<std::vector<int>>> mVarModFormula = ModificationUtils::getInstance().ExtVarModFormula(oligo.sSequence, vVarModEntrys);

        for (const auto& formulaMap : mVarModFormula) {
            if (this->mVarModIndex.find(formulaMap.first) != this->mVarModIndex.end()) {
                oligo.vVarModFormulaIndices.push_back(this->mVarModIndex[formulaMap.first]);
                oligo.mVarModPosition[this->mVarModIndex[formulaMap.first]] = formulaMap.second;
            }
            else {
                VarModFormula varModFormula;
                varModFormula.sFormula = formulaMap.first;
                int iIndex = vVarMods.size();
                this->mVarModIndex[formulaMap.first] = iIndex;
                vVarMods.push_back(varModFormula);
                oligo.vVarModFormulaIndices.push_back(iIndex);
                oligo.mVarModPosition[iIndex] = formulaMap.second;
            }
        }
    }

    // form detail information of varMods
    for (int i = 0; i < vVarMods.size(); i++)
    {
        VarModFormula& varMod = vVarMods[i];
        
        // while (pSearchThreadPool->jobs_.size() >= CommonValues::iMaxThreadProc)
        // {
        //     pSearchThreadPool->wait_for_available_thread();
        // }
        // pSearchThreadPool->doJob(std::bind(ModificationUtils::FormVarMod, std::ref(varMod), vModificationList));

        this->FormVarMod(varMod, vVarModEntrys);

        if (!(i % 1000))
        {
            WriteOutUtils::getInstance().ShowProgressBar(i, vVarMods.size(), "Form Variable Modification", "");
        }
    }

    // pSearchThreadPool->wait_on_threads();
    WriteOutUtils::getInstance().ShowProgressBar(vVarMods.size(), vVarMods.size(), "Form Variable Modification", "done\n");

    return true;

}

bool ModificationUtils::FormVarMod(VarModFormula& varMod, vector<VarModEntry> vVarModEntrys)
{   
    string sFormula = varMod.sFormula;                       //even sSequence == "", it should also work
    if (sFormula.empty())
    {
        varMod.vvVarModList.push_back(vector<Modification>());
        varMod.vVarModMassList.push_back(0.0);
        varMod.varModMassRange = DoubleRange(0.0, 0.0);
        return true;
    }

    // get all possible modification index for each nucleotide, for generate cartesian product
    vector<vector<int>> vvModIndexSet;
    for (const auto& nucleotide : sFormula)
    {
        vector<int> vModIndexSet;
        for (size_t iIndex = 0; iIndex < vVarModEntrys.size(); iIndex++)
        {
            const auto& varModEntry = vVarModEntrys[iIndex];
            if (varModEntry.sAllowNT.find(nucleotide) != string::npos)
            {
                vModIndexSet.push_back(iIndex);
            }
        }
        vvModIndexSet.push_back(vModIndexSet);
    }

    // calculate cartesian product of all possible modification index
    vector<vector<int>> vvCombinationSet;
    vector<int> currentCombination(sFormula.size());
    cartesianCombinations(vvModIndexSet, currentCombination, vvCombinationSet, 0);

    // check if the combination is legal
    for (auto& vCombination : vvCombinationSet)
    {
        std::map<int, int> mModEntryCount;
        bool bLegal = true;
        for (auto& index : vCombination)
        {
            mModEntryCount[index]++;
        }

        for (auto& entry : mModEntryCount)
        {
            if (entry.second > vVarModEntrys[entry.first].iMaxNum)
            {
                bLegal = false;
                break;
            }
        }
        if (!bLegal)
            continue;

        // store modification && calculate mass for each combination
        vector<Modification> vModificationCombination;
        double dMassSum = 0.0;
        for (auto& index : vCombination)
        {
            vModificationCombination.push_back(vVarModEntrys[index].modification);
            dMassSum += vVarModEntrys[index].modification.dModMass;
        }

        varMod.vvVarModList.push_back(vModificationCombination);
        varMod.vVarModMassList.push_back(dMassSum);
    }
    varMod.CalcMassRange();

    return true;
}

std::map<std::string, std::vector<std::vector<int>>> ModificationUtils::ExtVarModFormula(const string& sSequence, const vector<VarModEntry>& vVarModEntrys)
{
    std::map<std::string, std::vector<std::vector<int>>> vVarModFormula;
    vVarModFormula[""] = vector<vector<int>>{vector<int>()}; // no modification case

    if (sSequence.empty()) {
        return vVarModFormula;
    }

    std::string sSubSequence;
    std::vector<int> vSubSequenceIndex;
    for (size_t i = 0; i < sSequence.size(); i++) {
        for (const auto& varModEntry : vVarModEntrys) {
            if (varModEntry.sAllowNT.find(sSequence[i]) != std::string::npos) {
                sSubSequence += sSequence[i];
                vSubSequenceIndex.push_back(i);
                break; // Exit the loop once we find a match
            }
        }
    }

    std::map<char, int> mNucleotideCountMap; //{nucleotide, max_num}
    mNucleotideCountMap['A'] = mNucleotideCountMap['C'] = mNucleotideCountMap['G'] = mNucleotideCountMap['T'] = mNucleotideCountMap['U'] = 0;
    for (const auto& varModEntry : vVarModEntrys) {
        for (const auto& nucleotide : varModEntry.sAllowNT) {
            if (mNucleotideCountMap[nucleotide] < varModEntry.iMaxNum) {
                mNucleotideCountMap[nucleotide] = varModEntry.iMaxNum;
            }
        }
    }

    int iMaxVarModNum = min(Param::getInstance().options.iMaxVarMods, sSubSequence.size());
    // calculate all possible combination sets
    for (int iSelectNum = 1; iSelectNum <= iMaxVarModNum; iSelectNum++) {
        int iNCK = CombinatoricsUtils::getInstance().nChooseK(sSubSequence.size(), iSelectNum);
        int** vvCombinationSets = CombinatoricsUtils::getInstance().makeCombinations(sSubSequence.size(), iSelectNum, iNCK);
        for (int iCombinationIndex = 0; iCombinationIndex < iNCK; iCombinationIndex++) {
            int* vCombination = vvCombinationSets[iCombinationIndex];
            std::string sCombination;
            std::vector<int> vCombinationIndex;
            for (int i = 0; i < iSelectNum; i++) {
                sCombination += sSubSequence[vCombination[i]];
                vCombinationIndex.push_back(vSubSequenceIndex[vCombination[i]]);
            }
            
            // check if the sCombination is legal
            if (sCombination.size() <= 0)
                continue;
            
            bool bLegal = true;
            // scombination 中每个核苷酸的个数不能超过mNucleotideCountMap中对应的值
            std::map<char, int> mCombinationCountMap;
            mCombinationCountMap['A'] = mCombinationCountMap['C'] = mCombinationCountMap['G'] = mCombinationCountMap['T'] = mCombinationCountMap['U'] = 0;
            for (const auto& nucleotide : sCombination) {
                mCombinationCountMap[nucleotide]++;
            }
            for (const auto& entry : mCombinationCountMap) {
                if (entry.second > mNucleotideCountMap[entry.first]) {
                    bLegal = false;
                    break;
                }
            }
            
            // scombination 中对于每种可变修饰，其修饰总数不能超过vVarModEntrys中对应的值
            for (const auto& varModEntry : vVarModEntrys) {
                int iCount = 0;
                for (const auto& cNT: sCombination) {
                    if (varModEntry.sAllowNT.find(cNT) != std::string::npos) {
                        iCount++;
                    }
                }
                if (iCount > varModEntry.iMaxNum) {
                    bLegal = false;
                    break;
                }
            }

            if (bLegal) {
                if (vVarModFormula.find(sCombination) == vVarModFormula.end()) {
                    vVarModFormula[sCombination] = vector<vector<int>>{vCombinationIndex};
                }
                else {
                    vVarModFormula[sCombination].push_back(vCombinationIndex);
                }
            }

        }

        // delete combination sets
        for (int iCombinationIndex = 0; iCombinationIndex < iNCK; iCombinationIndex++) {
            delete[] vvCombinationSets[iCombinationIndex];
        }
        delete[] vvCombinationSets;

    }

    return vVarModFormula;
}
