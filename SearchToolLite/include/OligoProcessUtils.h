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


#ifndef _OLIGOPROCESSUTILS_H_
#define _OLIGOPROCESSUTILS_H_

#include "Common.h"
#include "DataStructure.h"
#include "Param.h"
#include "ThreadPool.h"

#include <map>
#include <vector>
#include <string>

class OligoProcessUtils
{
public:
   // struct OligoModInfo
   // {
   //    int iVarModFormulaIndex;
   //    int iPositionIndex;
   //    int iVarModIndex;

   //    OligoModInfo(int iVarModFormulaIndex, int iPositionIndex, int iVarModIndex)
   //    {
   //       this->iVarModFormulaIndex = iVarModFormulaIndex;
   //       this->iPositionIndex = iPositionIndex;
   //       this->iVarModIndex = iVarModIndex;  
   //    }

   //    OligoModInfo& operator=(const OligoModInfo& other)
   //    {
   //       if (this != &other)
   //       {
   //          iVarModFormulaIndex = other.iVarModFormulaIndex;
   //          iPositionIndex = other.iPositionIndex;
   //          iVarModIndex = other.iVarModIndex;
   //       }
   //       return *this;
   //    }
   // };

	OligoProcessUtils();
   ~OligoProcessUtils();

   bool ReadFastaAndDigest(ThreadPool *tp, vector<NucleicAcid>& vNucleicAcid, vector<Oligonucleotide>& vOligo); 

   bool ReadFasta(const std::string& fastaFilename, vector<NucleicAcid>& vNucleicAcid);

   void Digest(const NucleicAcid& nucleicAcid, vector<Oligonucleotide>& vOligo);

   vector<OligoFragment> GenerateFragments(const Oligonucleotide& oligo, int iCharge, const OligonucleotideID& oligoID);

   vector<OligoFragment> GenerateNLFragments(const Oligonucleotide& oligo, int iCharge, const OligonucleotideID& oligoID);

   bool GenerateDecoyOligo(vector<Oligonucleotide>& vOligo, int iDecoyOption);

   static OligoProcessUtils& getInstance()
   {
       static OligoProcessUtils instance;
       return instance;
   }

private:

   bool ReadEnzymeFormula(const std::string& sEnzymeDatabase, const std::string& sEnzymeName);

   std::map<string, std::vector<int>>  mOligoIndex;
   std::vector<std::string> vEnzymeFormula;
};

#endif
