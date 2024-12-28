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


#ifndef _CALCULATEUTILS_H_
#define _CALCULATEUTILS_H_

#include "Common.h"
#include "DataStructure.h"
#include "Param.h"

class CalculateUtils
{
public:
    CalculateUtils();
    ~CalculateUtils();

    static bool FastXcorrCalc(vector<double> theoreticalSignal, vector<double> experimentSignal, int offset, double& xcorr);

    static bool VectorDot(vector<double> vec1, vector<double> vec2, double& result);

    static int BinarySearchOligonucleotideByMass(const vector<Oligonucleotide>& vOligonucleotideList, double dMass);          //return the first oligonucleotide with mass larger than dMass

    static int BinarySearchOligonucleotideByLowestMass(const vector<Oligonucleotide>& vOligonucleotideList, double dMass);    //return the first oligonucleotide with lowest_mass larger than dMass

    static int BinarySearchOligonucleotideByHighestMass(const vector<Oligonucleotide>& vOligonucleotideList, double dMass);    //return the final oligonucleotide with uppest_mass smaller than dMass

    static void SortOligonucleotideByMass(vector<Oligonucleotide>& vOligonucleotideList);

    static void SortOligonucleotideByLowestMass(vector<Oligonucleotide>& vOligonucleotideList);

    static void SortOligonucleotideByHighestMass(vector<Oligonucleotide>& vOligonucleotideList);

    static void SortQueryResultBySp(Query& query);

    static void SortQueryByPrecursorMZ(vector<Query>& vQueryList);

private:

    static bool PreprocessSpectrumSignal(vector<double>& inputSignal, int offset);

};

class CombinatoricsUtils
{
public:
    CombinatoricsUtils();
    ~CombinatoricsUtils();

    static int** makeCombinations(int n, int r, int count);

    static int nChooseK(int n, int k);

    static int getCombinationCount(int n, int k);

    static void initBinomialCoefficients(const int n, const int k);

};

#endif // _COMETCALCULATEUTILS_H_
