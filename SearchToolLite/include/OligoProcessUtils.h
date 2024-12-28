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


#ifndef _SEQPROCESSUTILS_H_
#define _SEQPROCESSUTILS_H_

#include "Common.h"
#include "DataStructure.h"
#include "Param.h"
#include "ThreadPool.h"

class OligoProcessUtils
{
public:
	OligoProcessUtils();
   ~OligoProcessUtils();

   static bool ReadNucleicAcidAndDigest(ThreadPool *tp, vector<NucleicAcid>& vNucleicAcid, vector<Oligonucleotide>& vOligonucleotide); //parameter use global static parameter

   static bool ReadNucleicAcid(DBInfo& databaseInfo, vector<NucleicAcid>& vNucleicAcid);

   static void DigestNucleicAcid(NucleicAcid nucleicAcid, vector<Oligonucleotide>& vOligonucleotide);

   static vector<OligonucleotideFragment> GenerateOligonucleotideFragment(const Oligonucleotide& oligonucleotide, IonSeries::Ions ionType, int iCharge, vector<pair<int, Modification>> vVarModList);

   static vector<OligonucleotideFragment> GenerateNLOligonucleotideFragment(const Oligonucleotide& oligonucleotide, int iCharge);

   static bool GenerateDecoyOligonucleotide(vector<Oligonucleotide>& vOligonucleotideList, int iDecoyOption);

private:

   static bool CheckEnzymeTermini(string seq, int iStart, int iEnd, Enzyme& enzyme);

   static unordered_map<string, int>  decoyDuplicateMap;
};

#endif
