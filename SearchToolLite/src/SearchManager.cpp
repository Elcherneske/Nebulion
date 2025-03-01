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

#include "SearchManager.h"
#include "MassUtils.h"
#include "OligoProcessUtils.h"
#include "SpectrumPreprocess.h"
#include "ModificationUtils.h"
#include "WriteOutUtils.h"
#include "CalculateUtils.h"
#include "SearchUtils.h"

/******************************************************************************
*
* Static helper functions
*
******************************************************************************/

static InputFile::InputType GetInputType(string sFileName)
{
   int iLen = sFileName.size();

   if (sFileName.substr(iLen - 6, 6) == ".mzXML"
      || sFileName.substr(iLen - 5, 5) == ".mzML"
      || sFileName.substr(iLen - 9, 9) == ".mzXML.gz"
      || sFileName.substr(iLen - 8, 8) == ".mzML.gz")
      return InputFile::InputType_MZXML;
   else if (sFileName.substr(iLen - 4, 4) == ".raw")
   {
      return InputFile::InputType_RAW;
   }
   else if (sFileName.substr(iLen - 4, 4) == ".ms2"
         || sFileName.substr(iLen - 5, 5) == ".cms2")
   {
      return InputFile::InputType_MS2;
   }
   else if (sFileName.substr(iLen - 4, 4) == ".mgf")
   {
      return InputFile::InputType_MGF;
   }

   return InputFile::InputType_UNKNOWN;
}


static void SetMSLevelFilter(MSReader &mstReader)
{
   vector<MSSpectrumType> msLevel;

   if (Param::getInstance().options.iMSLevel == 3)
      msLevel.push_back(MS3);
   else
      msLevel.push_back(MS2);

   mstReader.setFilter(msLevel);
}

/******************************************************************************
*
* SearchManager class implementation.
*
******************************************************************************/

SearchManager::SearchManager()
{
   NebulionInterface::_tp = new ThreadPool();
}

SearchManager::~SearchManager()
{
   if (_tp != NULL)
      delete _tp;
   _tp = NULL;
}

bool SearchManager::IsSearchError()
{
    //return Param::g_cometStatus.IsError();
    return false;
}

void SearchManager::GetStatusMessage(string &strStatusMsg)
{
   //Param::g_cometStatus.GetStatusMsg(strStatusMsg);
}

void SearchManager::CancelSearch()
{
   //Param::g_cometStatus.SetStatus(Result_Cancelled, string("Search was cancelled."));
}

bool SearchManager::IsCancelSearch()
{
   //return Param::g_cometStatus.IsCancel();
    return false;
}

void SearchManager::ResetSearchStatus()
{
   //Param::g_cometStatus.ResetStatus();
}

bool SearchManager::DoSearch()
{
    bool bSucceeded = true;
    ThreadPool *tp = NebulionInterface::_tp;
    tp->fillPool(Param::getInstance().options.iThreadNum <= 0 ? 0 : Param::getInstance().options.iThreadNum - 1);

   // read fasta and digest to oligonucleotide
   bSucceeded = OligoProcessUtils::getInstance().ReadFastaAndDigest(tp, Param::getInstance().vNucleicAcids, Param::getInstance().vOligos);
   if (!bSucceeded)
      return bSucceeded;

   // add decoy oligonucleotide
   if (Param::getInstance().options.bDecoySearch)  //judge whether to add decoy
   {
       cout << "Begin to add decoy Oligonucleotides" << endl << endl;
       bSucceeded = OligoProcessUtils::getInstance().GenerateDecoyOligo(Param::getInstance().vOligos, Param::getInstance().options.iDecoyKind);
       if (!bSucceeded)
           return bSucceeded;
   }

   // add various modification
   cout << "begin add various modification" << endl << endl;
   bSucceeded = ModificationUtils::getInstance().AddVarMod(_tp, Param::getInstance().vOligos, Param::getInstance().vVarModFormulas, Param::getInstance().options.sVarModsDatabase);
   if (!bSucceeded)
      return bSucceeded;

   WriteOutUtils::getInstance().OligosCSVOutput("D:/code/VSC++/Nebulion/resource/oligos.csv", Param::getInstance().vOligos);

   // bin oligonucleotide by mass
   cout << "begin bin oligonucleotide by mass" << endl << endl;
   bSucceeded = PrecursorMassIndexUtils::getInstance().GeneratePrecursorMassIndex(Param::getInstance().vOligos);
   if (!bSucceeded)
      return bSucceeded;

   for (size_t iFileIndex = 0; iFileIndex < Param::getInstance().vInputFiles.size(); ++iFileIndex) //对每个inputfile读取并搜索
   {
      InputFile& inputFile = Param::getInstance().vInputFiles[iFileIndex];
      vector<ExpSpectrum> vSpectrums;
      vector<Query> vQueries;

      cout << " - Input file: " + string(inputFile.sFileName) << endl;
      // load spectrum && preprocess
      SetMSLevelFilter(SpectrumPreprocessUtils::getInstance().GetMSReader());
      bSucceeded = SpectrumPreprocessUtils::getInstance().LoadSpectrum(inputFile, vSpectrums);

      // limit by precursor mz match && form query
      for (size_t iSpecIndex = 0; iSpecIndex < vSpectrums.size(); iSpecIndex++)
      {
          ExpSpectrum& expspectrum = vSpectrums[iSpecIndex];
          Query query(expspectrum);

          bSucceeded = SearchUtils::getInstance().PrecursorFilter(tp, query, Param::getInstance().vOligos);
          bSucceeded = SearchUtils::getInstance().RoughScoring(tp, query);

          if (!bSucceeded)
          {
              std::cout << "error" << std::endl;
              continue;
          }

          if (query.vOligoEntries.size() > 0)
          {
            vQueries.push_back(query);
          }

          WriteOutUtils::getInstance().ShowProgressBar(iSpecIndex, vSpectrums.size(), "Processing spectrum", " ");
      }
      
      WriteOutUtils::getInstance().ShowProgressBar(vSpectrums.size(), vSpectrums.size(), "Processing spectrum", "done\n");
      {
          cout << "MS1 query match :" << vQueries.size() << endl;
      }
      
      WriteOutUtils::getInstance().QueriesCSVOutput("D:/code/VSC++/Nebulion/resource/result.csv", vQueries);

      cout << " - done." << endl << endl;
      
   }

   return true;
}
