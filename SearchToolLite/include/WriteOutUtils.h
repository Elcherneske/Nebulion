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

#ifndef _WRITEOUTUTILS_
#define _WRITEOUTUTILS_

#include "Common.h"
#include "DataStructure.h"

class WriteOutUtils
{
public:
	WriteOutUtils();
   ~WriteOutUtils();

   static bool SaveQueries(string outputFilename, vector<Query>& vQueryList);

   static bool OutputParams(ofstream& outputStream);

   static void ShowProgressBar(double iProgress, double iTotal, string sPreDescribe, string sPostDescribe);

   static bool QueryCSVOutput(ofstream& outputStream, Query& query, bool head);

   static bool TmpQueryOut(ofstream& outputStream, Query& query);

   static void OutputOligonucleotides(ofstream& outputStream, vector<Oligonucleotide>& vOligonucleotideList);

   static bool TmpSpectrumOut(ofstream& outputStream, vector<ExpSpectrum>& vExpSpectrumList);

private:

};

#endif
