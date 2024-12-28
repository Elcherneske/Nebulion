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


#ifndef _SEARCHUTILS_H_
#define _SEARCHUTILS_H_

#include "Common.h"
#include "DataStructure.h"
#include "Param.h"

class SearchUtils
{
public:
	SearchUtils();
   ~SearchUtils();

   static bool MS1Search(Query& query, vector<Oligonucleotide>& vOligonucleotideList);
   static bool MS1Search(ThreadPool* _tp, vector<Query>& vQueries, const vector<Oligonucleotide>& vOligonucleotideList);

   static bool QueryScoring(ThreadPool* _tp, Query& query, const vector<Oligonucleotide>& vOligonucleotide);
   static bool QueryScoring(ThreadPool* _tp, vector<Query>& vQueries, const vector<Oligonucleotide>& vOligonucleotide);

   static bool EntryScoring(const CometSpectrum& spec, const Oligonucleotide& oligonucleotide, const int iWhichVarModCombination, Scores& result);

private:
   // Core search functions
   static bool GenerateTheoreticalSpectrum(vector<OligonucleotideFragment> fragmentList, int iArraySize, int iMaxFragmentCharge, vector<double>& spectrum);

   static bool GenerateExperimentalSpectrum(const CometSpectrum& cometSpectrum, vector<double>& spectrum);

   static vector<pair<int, int>> FragmentMatchCount(const CometSpectrum& cometSpectrum, const vector<OligonucleotideFragment>& vFragmentList);

};

#endif
