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
// static void GetHostName()
// {
// #ifdef _WIN32
//    WSADATA WSAData;
//    WSAStartup(MAKEWORD(1, 0), &WSAData);

//    if (gethostname(g_staticParams.szHostName, CommonValues::iFileBufferSize) != 0)
//       strcpy(g_staticParams.szHostName, "locahost");

//    WSACleanup();
// #else
//    if (gethostname(g_staticParams.szHostName, CommonValues::iFileBufferSize) != 0)
//       strcpy(g_staticParams.szHostName, "locahost");
// #endif

//    char *pStr;
//    if ((pStr = strchr(g_staticParams.szHostName, '.'))!=NULL)
//       *pStr = '\0';
// }

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

   if (Param::g_staticParams.options.iMSLevel == 3)
      msLevel.push_back(MS3);
   else
      msLevel.push_back(MS2);

   mstReader.setFilter(msLevel);
}



static bool ValidateSequenceDatabaseFile()
{
   ifstream fsCheck;

   // if .idx database specified but does not exist, first see if corresponding
   // fasta exists and if it does, create the .idx file
   if (Param::g_staticParams.databaseInfo.sDatabase.find(".idx") != string::npos)
   {
      fsCheck = ifstream(Param::g_staticParams.databaseInfo.sDatabase);
      if (!fsCheck) //not .idx file
      {
         fsCheck = ifstream(Param::g_staticParams.databaseInfo.sFilename);
         if (!fsCheck)
         {
            string strErrorMsg = " Error - peptide index file \"" + Param::g_staticParams.databaseInfo.sDatabase + "\" and corresponding FASTA file\n         are both missing.\n";
            logerr(strErrorMsg.c_str());
            return false;
         }
         else
            return true;
      }
      else
      {
         fsCheck = ifstream(Param::g_staticParams.databaseInfo.sFilename);
         if (!fsCheck)
         {
            string strErrorMsg = " Error - peptide index file \"" + Param::g_staticParams.databaseInfo.sDatabase + "\" specified is present\n   but corresponding FASTA file\n         are both missing.\n";
            logerr(strErrorMsg.c_str());
            return false;
         }
         else
            return true;
      }
   }
   else
   {
      fsCheck = ifstream(Param::g_staticParams.databaseInfo.sDatabase);
      if (!fsCheck)
      {
         string strErrorMsg = " Error - Database file \"" + Param::g_staticParams.databaseInfo.sDatabase + "\" is missing.\n";
         logerr(strErrorMsg.c_str());
         return false;
      }
      else
      {
         return true;
      }
   }

   return true;
}

static bool ValidateScanRange()
{
   if (Param::g_staticParams.options.scanRange.iEnd < Param::g_staticParams.options.scanRange.iStart && Param::g_staticParams.options.scanRange.iEnd != 0)
   {
      char szErrorMsg[512];
      sprintf(szErrorMsg, " Error - start scan is %d but end scan is %d.\n The end scan must be >= to the start scan.\n",
            Param::g_staticParams.options.scanRange.iStart,
            Param::g_staticParams.options.scanRange.iEnd);
      string strErrorMsg(szErrorMsg);
      // g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
      logerr(szErrorMsg);
      return false;
   }

   return true;
}

static bool ValidatePeptideLengthRange()
{
   if (Param::g_staticParams.options.oligonucleotideLengthRange.iEnd < Param::g_staticParams.options.oligonucleotideLengthRange.iStart 
      && Param::g_staticParams.options.oligonucleotideLengthRange.iEnd != 0)
   {
      char szErrorMsg[512];
      sprintf(szErrorMsg, " Error - peptide length range set as %d to %d.\n The maximum length must be >= to the minimum length.\n",
            Param::g_staticParams.options.oligonucleotideLengthRange.iStart,
            Param::g_staticParams.options.oligonucleotideLengthRange.iEnd);
      string strErrorMsg(szErrorMsg);
      // g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
      logerr(szErrorMsg);
      return false;
   }

   return true;
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

void SearchManager::SetOutputFileBaseName(string pszBaseName)
{

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
    tp->fillPool(Param::g_staticParams.options.iNumThreads <= 0 ? 0 : Param::g_staticParams.options.iNumThreads - 1);
    Param::initMutex();

    if (!Param::CheckParamsValidation()) return false;

    {
        cout << "Comet version \"" + Param::g_sCometVersion + "\"" << endl << endl;
    }

   {
        cout << "Begin read nuclic acid and digest to oligoncleotides" << endl << endl;
   }

   bSucceeded = OligoProcessUtils::ReadNucleicAcidAndDigest(tp, Param::g_vNucleicAcid, Param::g_vOligonucleotide); 
   if (!bSucceeded)
      return bSucceeded;

   if (Param::g_staticParams.options.iDecoySearch)  //judge whether to add decoy
   {
       cout << "Begin to add decoy Oligonucleotides" << endl << endl;
       bSucceeded = OligoProcessUtils::GenerateDecoyOligonucleotide(Param::g_vOligonucleotide, Param::g_staticParams.options.iDecoyKind);
       if (!bSucceeded)
           return bSucceeded;
   }
   
   {
       ofstream seqOutputStream("sequence.txt");
       WriteOutUtils::OutputOligonucleotides(seqOutputStream, Param::g_vOligonucleotide);
   }

   {
       cout << "begin add various modification" << endl << endl;
   }
   bSucceeded = ModificationUtils::AddVarMod(_tp, Param::g_vOligonucleotide, Param::g_staticParams.modificationInfo.sVariableModificationFilename);
   if (!bSucceeded)
      return bSucceeded;

   CalculateUtils::SortOligonucleotideByHighestMass(Param::g_vOligonucleotide);

   for (int i=0; i<(int)Param::g_vInputFiles.size(); ++i) //对每个inputfile读取并搜索
   {
      InputFile inputFile = Param::g_vInputFiles[i];
      vector<CometSpectrum> vSpectrumList;
      vector<Query> vQueries;

      {
          cout << " - Input file: " + string(inputFile.sFileName) << endl;
          SetMSLevelFilter(SpectrumPreprocessUtils::getMSReader());
          // Load and preprocess all the spectra.
          {
              cout << "\n   - Load spectra and preprocess:" << endl;
          }
          bSucceeded = SpectrumPreprocessUtils::LoadSpectrum(inputFile, vSpectrumList);
      }

      for (int spec_index = 0; spec_index < vSpectrumList.size(); spec_index++)
      {
          CometSpectrum& cometspectrum = vSpectrumList[spec_index];
          Query query;
          query.cometSpectrum = &cometspectrum;

          bSucceeded = SearchUtils::MS1Search(query, Param::g_vOligonucleotide);

          if (!bSucceeded)
          {
              cout << "error" << endl;
              continue;
          }

          if (query.vResults.size() != 0)
          {
              vQueries.push_back(query);
          }
      }

      {
          cout << "MS1 query match :" << vQueries.size() << endl;
      }

      for (int query_index = 0; query_index < vQueries.size(); query_index++)
      {
          Query& query = vQueries[query_index];
          bSucceeded = SearchUtils::QueryScoring(_tp, query, Param::g_vOligonucleotide);
          CalculateUtils::SortQueryResultBySp(query);

          WriteOutUtils::ShowProgressBar(query_index, vQueries.size(), "MS2 Scoring:", "");
      }
      WriteOutUtils::ShowProgressBar(vQueries.size(), vQueries.size(), "MS2 Scoring:", "done \n");

      CalculateUtils::SortQueryByPrecursorMZ(vQueries);
         
      int iFirst = 0;
      ofstream csvFile("result.csv");
      for (Query& query : vQueries)
      {
          if (iFirst == 0)
          {
              WriteOutUtils::QueryCSVOutput(csvFile, query, true);
              iFirst++;
          }
          else
              WriteOutUtils::QueryCSVOutput(csvFile, query, false);
      }

      ofstream spectrumOutputStream("spectrum.txt");
      WriteOutUtils::TmpSpectrumOut(spectrumOutputStream, vSpectrumList);

      cout << " - done." << endl << endl;
      
   }

   return true;
}
