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


#include "WriteOutUtils.h"
#include "Common.h"
#include "DataStructure.h"
#include "Param.h"
#include "MassUtils.h"

#include <fstream>
#include <string>
#include <vector>


WriteOutUtils::WriteOutUtils()
{
}


WriteOutUtils::~WriteOutUtils()
{
}

void WriteOutUtils::ShowProgressBar(double iProgress, double iTotal, const std::string& sPreDescribe, const std::string& sPostDescribe) 
{
    // 计算进度百分比
    int percentage = (int)((iProgress * 100) / iTotal);

    // 进度条的长度
    const int barWidth = 50;
    int pos = barWidth * iProgress / iTotal;

    // 在同一行显示进度条，\r 会使输出回到行首，从而覆盖之前的内容
    std::cout << "\r" + sPreDescribe + "[";
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) {
            std::cout << "=";  // 已完成的部分
        }
        else if (i == pos) {
            std::cout << ">";  // 当前进度
        }
        else {
            std::cout << " ";  // 未完成的部分
        }
    }
    std::cout << "] " << percentage << "%" << "\t" << sPostDescribe;
    std::cout.flush();  // 刷新输出，确保进度条实时更新
}

bool WriteOutUtils::SaveQueries(const std::string& outputFilename, const std::vector<Query>& vQueryList)
{
   ofstream outputStream(outputFilename);

   for(auto& query : vQueryList)
   {
      // WriteOutUtils::WriteSingleQuery(outputStream, query);
   }

   return true;
}

void WriteOutUtils::QueriesCSVOutput(const std::string& outputFilename, const std::vector<Query>& vQueries, int iTopN)
{
   ofstream outputStream(outputFilename);

   if(!outputStream)
   {
      std::cerr << "Error: wrong outputfile stream" << std::endl;
      return;
   }

   outputStream << "Precursor Mz, Sequence, Mass, NucleicAcid Origin, IsDecoy, 3-end modificastion, 5-end modification, Static Modifications, Variable Modification, XCorr\n";

   for (const auto& query : vQueries) 
   {
      int iEntryNum = min(iTopN, (int)query.vOligoEntries.size());

      for (int iEntryIndex = 0; iEntryIndex < iEntryNum; iEntryIndex++)
      {
         const auto& entry = query.vOligoEntries[iEntryIndex];
         const Oligonucleotide& oligo = Param::getInstance().vOligos[entry.oligoID.iOligoIndex];
         outputStream << fixed << setprecision(8);
        
         // 输出 Precursor mz
         outputStream << query.expSpectrum.dExpMZ << ",";
         
         // 输出 Oligonucleotide 信息
         outputStream << oligo.sSequence << ",";
         outputStream << MassUtils::getInstance().CalOligoMass(oligo, entry.oligoID, true, true, true) << ",";
         std::string sNAPos;
         for (const auto& nucleicAcidIndex : oligo.vWhichNucleicAcid)
         {
            if (nucleicAcidIndex.iWhichNucleicAcid >= 0 && nucleicAcidIndex.iWhichNucleicAcid < Param::getInstance().vNucleicAcids.size())
            {
               sNAPos += Param::getInstance().vNucleicAcids[nucleicAcidIndex.iWhichNucleicAcid].sName + ":";
               sNAPos += to_string(nucleicAcidIndex.iNABeginPos)+ "-" + to_string(nucleicAcidIndex.iNAEndPos) + ";";
            }
            else
            {
               sNAPos += "-;";
               break;
            }
         }
         outputStream << sNAPos << ",";
         outputStream << (oligo.bDecoy ? "Yes" : "No") << ",";

         outputStream << oligo.end3TermMod.toString() << ",";
         outputStream << oligo.end5TermMod.toString() << ",";
         
         // 输出 Static Modifications
         std::string sModsStr;
         for (const auto& mod : oligo.vMods) {
            sModsStr += to_string(mod.first) + ":" + mod.second.toString() + "; "; 
         }
         outputStream << sModsStr << ",";

         // 输出variable modification
         outputStream << "-" << ",";

         // 输出打分信息
         outputStream << entry.scores.xCorr << "\n";
      }
   
   }
}

void WriteOutUtils::OligosCSVOutput(const std::string& outputFilename, const std::vector<Oligonucleotide>& vOligos)
{  
   ofstream outputStream(outputFilename);

   if(!outputStream)
   {
      std::cerr << "Error: wrong outputfile stream" << std::endl;
      return;
   }

   outputStream << std::setprecision(8);

   outputStream << "Sequence, Nuclic Acid, Decoy, 3-end modification, 5-end modification, Static Modifications, Variable Modifications" << endl;

   // 输出每个 oligonucleotide 的信息
   for (const auto& oligo : vOligos)
   {
       outputStream << oligo.sSequence << ",";
       std::string sNAPos;
      for (const auto& nucleicAcidIndex : oligo.vWhichNucleicAcid)
      {
         if (nucleicAcidIndex.iWhichNucleicAcid >= 0 && nucleicAcidIndex.iWhichNucleicAcid < Param::getInstance().vNucleicAcids.size())
         {
            sNAPos += Param::getInstance().vNucleicAcids[nucleicAcidIndex.iWhichNucleicAcid].sName + ":";
            sNAPos += to_string(nucleicAcidIndex.iNABeginPos)+ "-" + to_string(nucleicAcidIndex.iNAEndPos) + ";";
         }
         else
         {
            sNAPos += "-;";
            break;
         }
      }
       outputStream << sNAPos << ",";
       outputStream << (oligo.bDecoy ? "Yes" : "No") << ",";

       outputStream << oligo.end3TermMod.toString() << ",";
       outputStream << oligo.end5TermMod.toString() << ",";

       // 输出 Static Modifications
       std::string sModsStr;
       for (const auto& mod : oligo.vMods) {
           sModsStr += to_string(mod.first) + ":" + mod.second.toString() + "; ";
       }
       outputStream << sModsStr << ",";

       // 输出variable modification
       outputStream << "-" << endl;
   }

}

