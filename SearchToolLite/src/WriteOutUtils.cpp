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


WriteOutUtils::WriteOutUtils()
{
}


WriteOutUtils::~WriteOutUtils()
{
}

void WriteOutUtils::ShowProgressBar(double iProgress, double iTotal, string sPreDescribe, string sPostDescribe) 
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

bool WriteOutUtils::SaveQueries(string outputFilename, vector<Query>& vQueryList)
{
   ofstream outputStream(outputFilename);

   WriteOutUtils::OutputParams(outputStream);

   for(auto& query : vQueryList)
   {
      // WriteOutUtils::WriteSingleQuery(outputStream, query);
   }

   return true;
}

bool WriteOutUtils::QueryCSVOutput(ofstream& outputStream, Query& query, bool head)
{
   if(!outputStream)
      return false;
   
   if (head)
      outputStream << "Precursor Mz,Sequence,Mass,NucleicAcid Origin,IsDecoy,3-end modificastion,5-end modification,Static Modifications,Variable Modification,XCorr,Sp,Cn,Expect,MatchedIons,TotalIons\n";

   for (const auto& result : query.vResults) 
   {
      if (query.cometSpectrum == nullptr)
      {
         string sErrorMsg = "Error: query cometSpectrum pointer is NULL; \n";
         logerr(sErrorMsg.c_str());
         continue;
      }

      const Oligonucleotide& oligonucleotide = result.oligo;

      double precursorMz = 0.0; 
      if (query.cometSpectrum) {
         precursorMz = query.cometSpectrum->spectrum.getMZ();
      }

      outputStream << fixed << setprecision(8);
      
      // 输出 Precursor mz
      outputStream << precursorMz << ",";
      
      // 输出 Oligonucleotide 信息
      outputStream << oligonucleotide.sSequence << ",";
      outputStream << oligonucleotide.dMass << ",";
      string NAPos;
      for (auto& info : oligonucleotide.whichNucleicAcid)
      {
         if (info.iWhichNucleicAcid >= 0 && info.iWhichNucleicAcid < Param::g_vNucleicAcid.size())
         {
            NAPos += Param::g_vNucleicAcid[info.iWhichNucleicAcid].sName + ":";
            NAPos += to_string(info.iNABeginPos)+ "-" + to_string(info.iNAEndPos) + ";";
         }
         else
         {
            NAPos = "-";
            break;
         }
      }
      outputStream << NAPos << ",";
      outputStream << (oligonucleotide.bDecoy ? "Yes" : "No") << ",";

      outputStream << oligonucleotide.end3TermMod.toString() << ",";
      outputStream << oligonucleotide.end5TermMod.toString() << ",";
      
      // 输出 Static Modifications
      string modsStr;
      for (const auto& mod : oligonucleotide.vMods) {
         modsStr += to_string(mod.first) + ":" + mod.second.toString() + "; "; 
      }
      outputStream << modsStr << ",";

      // 输出variable modification
      outputStream << "-" << ",";

      // 输出打分信息
      outputStream << result.scores.xCorr << ",";
      outputStream << result.scores.dSp << ",";
      outputStream << result.scores.dCn << ",";
      outputStream << result.scores.dExpect << ",";
      outputStream << result.scores.matchedIons << ",";
      outputStream << result.scores.totalIons << "\n";
   }
}

bool WriteOutUtils::OutputParams(ofstream& outputStream)
{
   if(!outputStream)
   {  
      string sErrorMsg = "Error: wrong outputfile stream";
      logerr(sErrorMsg.c_str());
      return false;
   }
   
   return true;
}

void WriteOutUtils::OutputOligonucleotides(ofstream& outputStream, vector<Oligonucleotide>& vOligonucleotideList)
{
   std::sort(vOligonucleotideList.begin(), vOligonucleotideList.end(), 
            [](const Oligonucleotide& a, const Oligonucleotide& b) 
            {return a.dMass < b.dMass;}
   );
   
   outputStream << std::setprecision(8);

   outputStream << "Sequence,Mass,Nuclic Acid,Decoy,3TermMod,5TermMod,Static Modification,Variable Modification" << endl;

   // 输出每个 oligonucleotide 的信息
   for (auto& oligonucleotide : vOligonucleotideList)
   {
       outputStream << oligonucleotide.sSequence << ",";
       outputStream << oligonucleotide.dMass << ",";
       string NAPos;
       for (auto& info : oligonucleotide.whichNucleicAcid)
       {
           if (info.iWhichNucleicAcid >= 0 && info.iWhichNucleicAcid < Param::g_vNucleicAcid.size())
           {
               NAPos += Param::g_vNucleicAcid[info.iWhichNucleicAcid].sName + ":";
               NAPos += to_string(info.iNABeginPos) + "-" + to_string(info.iNAEndPos) + ";";
           }
           else
           {
               NAPos = "-";
               break;
           }
       }
       outputStream << NAPos << ",";
       outputStream << (oligonucleotide.bDecoy ? "Yes" : "No") << ",";

       outputStream << oligonucleotide.end3TermMod.toString() << ",";
       outputStream << oligonucleotide.end5TermMod.toString() << ",";

       // 输出 Static Modifications
       string modsStr;
       for (const auto& mod : oligonucleotide.vMods) {
           modsStr += to_string(mod.first) + ":" + mod.second.toString() + "; ";
       }
       outputStream << modsStr << ",";

       // 输出variable modification
       outputStream << "-" << endl;
   }

}

bool WriteOutUtils::TmpSpectrumOut(ofstream& outputStream, vector<CometSpectrum>& vCometSpectrumList)
{
   outputStream << std::setprecision(16);
   for(auto& cometSpectrum : vCometSpectrumList)
   {
      outputStream << cometSpectrum.spectrum.getMZ() << "\t" << endl;
   }
   return true;
}

