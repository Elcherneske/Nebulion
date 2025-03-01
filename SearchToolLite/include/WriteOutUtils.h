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

   bool SaveQueries(const std::string& outputFilename, const std::vector<Query>& vQueryList);

   void ShowProgressBar(double iProgress, double iTotal, const std::string& sPreDescribe, const std::string& sPostDescribe);

   void QueriesCSVOutput(const std::string& outputFilename, const std::vector<Query>& vQueries, int iTopN = 5);

   bool TmpQueryOut(ofstream& outputStream, Query& query);

   void OligosCSVOutput(const std::string& outputFilename, const std::vector<Oligonucleotide>& vOligos);

   static WriteOutUtils& getInstance()
   {
       static WriteOutUtils instance;
       return instance;
   }

private:

};

#endif
