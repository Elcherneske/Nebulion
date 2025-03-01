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

    double FastXcorrCalculate(const vector<double>& theoreticalSignal, const vector<double>& experimentSignal, int offset);
    double FastXcorrCalculate(int* theoreticalSignal, int* experimentalSignal, int iArraySize, int offset);

    double VectorDot(const vector<double>& vec1, const vector<double>& vec2);
    double VectorDot(int* vec1, int* vec2, int iArraySize);

    void SortQueryResultBySp(Query& query);

    void SortQueryByPrecursorMZ(vector<Query>& vQueryList);

    static CalculateUtils& getInstance()
    {
        static CalculateUtils instance;
        return instance;
    }

private:

    vector<double> PreprocessSpectrumSignal(const vector<double>& inputSignal, int offset);
    void PreprocessSpectrumSignal(int* inputSignal, int iArraySize, int offset);

};

class CombinatoricsUtils
{
public:
    int** makeCombinations(int n, int r, int count);

    int nChooseK(int n, int k);

    int getCombinationCount(int n, int k);

    void initBinomialCoefficients(const int n, const int k);

    static CombinatoricsUtils& getInstance()
    {
        static CombinatoricsUtils instance;
        return instance;
    }

};

#endif // _COMETCALCULATEUTILS_H_
