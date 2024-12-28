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

#include "CalculateUtils.h"


CalculateUtils::CalculateUtils()
{
    
}
CalculateUtils::~CalculateUtils()
{

}

bool CalculateUtils::PreprocessSpectrumSignal(vector<double>& inputSignal, int offset)
{
    // Make fast xcorr spectrum.
    double dSum=0.0;
    int iWindowRange = 2 * offset + 1;
    double dInv = 1.0 / (iWindowRange - 1.0);
    vector<double> rawSignal = inputSignal;

    if (offset <= 0)
        return false;
    
    if (rawSignal.size() < offset)    //inputSignal is small
    {
        dSum = accumulate(rawSignal.begin(), rawSignal.end(), dSum);
        for (auto& it : inputSignal)
            it = dSum * dInv;
    }
    else
    {
        for (int i=0; i < offset; ++i)
            dSum += rawSignal[i];
        for (int i = offset; i < rawSignal.size() + offset; ++i)
        {
            if (i < rawSignal.size())
                dSum += rawSignal[i];
            if (i >= iWindowRange)
                dSum -= rawSignal[i-iWindowRange];
            inputSignal[i-offset] = (dSum - rawSignal[i-offset]) * dInv;
        }
    }
    return true;
}

bool CalculateUtils::FastXcorrCalc(vector<double> theoreticalSignal, vector<double> experimentSignal, int offset, double& xcorr)
{
    int max_size = max(theoreticalSignal.size(), experimentSignal.size());
    if (theoreticalSignal.size() < max_size) {
        theoreticalSignal.resize(max_size, 0.0);
    }

    if (experimentSignal.size() < max_size) {
        experimentSignal.resize(max_size, 0.0);
    }

    if(!CalculateUtils::PreprocessSpectrumSignal(experimentSignal, offset))
        return false;
    
    if(!VectorDot(theoreticalSignal, experimentSignal, xcorr))
        return false;

    return true;
}

bool CalculateUtils::VectorDot(vector<double> vec1, vector<double> vec2, double& result)
{
    if (vec1.size() != vec2.size() || vec1.size() * vec2.size() == 0)
        return false;
    
    result = inner_product(vec1.begin(), vec1.end(), vec2.begin(), 0.0);
    return true;
}   

void CalculateUtils::SortOligonucleotideByMass(vector<Oligonucleotide>& vOligonucleotideList)
{
    std::sort(vOligonucleotideList.begin(), vOligonucleotideList.end(), 
            [](const Oligonucleotide& a, const Oligonucleotide& b) 
            {return a.dMass < b.dMass;}
    );
}

void CalculateUtils::SortOligonucleotideByLowestMass(vector<Oligonucleotide>& vOligonucleotideList)
{
    std::sort(vOligonucleotideList.begin(), vOligonucleotideList.end(), 
            [](const Oligonucleotide& a, const Oligonucleotide& b) 
            {return (a.dMass + (Param::g_vVarModification[a.iVarModificationIndex]).varModMassRange.dStart) < (b.dMass + (Param::g_vVarModification[b.iVarModificationIndex]).varModMassRange.dStart);}
    );
}

void CalculateUtils::SortOligonucleotideByHighestMass(vector<Oligonucleotide>& vOligonucleotideList)
{
    std::sort(vOligonucleotideList.begin(), vOligonucleotideList.end(),
        [](const Oligonucleotide& a, const Oligonucleotide& b)
        {return (a.dMass + (Param::g_vVarModification[a.iVarModificationIndex]).varModMassRange.dEnd) < (b.dMass + (Param::g_vVarModification[b.iVarModificationIndex]).varModMassRange.dEnd);}
    );
}

void CalculateUtils::SortQueryResultBySp(Query& query) //默认由大到小
{
    std::sort(query.vResults.begin(), query.vResults.end(), 
        [](const Query::QueryEntry& a, const Query::QueryEntry& b)
        {
            if (a.scores.dSp > b.scores.dSp)
                return true;
            else if (a.scores.dSp < b.scores.dSp)
                return false;
            else
            {
                if(a.scores.xCorr > b.scores.xCorr)
                    return true;
                else if (a.scores.xCorr < b.scores.xCorr)
                    return false;
                else
                {
                    if (a.scores.matchedIons > b.scores.matchedIons)
                        return true;
                    else if (a.scores.matchedIons < b.scores.matchedIons)
                        return false;
                    else
                        return false;
                }
            }
        }
    );
}

void CalculateUtils::SortQueryByPrecursorMZ(vector<Query>& vQueryList)
{
    std::sort(vQueryList.begin(), vQueryList.end(), 
            [](const Query& a, const Query& b) 
            {
                return a.cometSpectrum->spectrum.getMZ() < b.cometSpectrum->spectrum.getMZ();
            }
    );
}

int CalculateUtils::BinarySearchOligonucleotideByMass(const vector<Oligonucleotide>& vOligonucleotideList, double dMass)
{
    int left = 0;
    int right = vOligonucleotideList.size() - 1;

    if (right < 0) {
        return -1;
    }

    while (left <= right) {
        int mid = left + (right - left) / 2;
    
        if (vOligonucleotideList[mid].dMass > dMass) 
            right = mid - 1;
        else 
            left = mid + 1;
        
    }

    if (left < vOligonucleotideList.size() && vOligonucleotideList[left].dMass > dMass) {
        return left;
    }

    return -1;
}

int CalculateUtils::BinarySearchOligonucleotideByLowestMass(const vector<Oligonucleotide>& vOligonucleotideList, double dMass) //should sort by lowest mass first
{
    int left = 0;
    int right = vOligonucleotideList.size() - 1;

    if (right < 0) 
    {
        return -1;
    }

    while (left <= right) 
    {
        int mid = left + (right - left) / 2;
    
        if (vOligonucleotideList[mid].dMass + (Param::g_vVarModification[vOligonucleotideList[mid].iVarModificationIndex]).varModMassRange.dStart > dMass) 
            right = mid - 1;
        else 
            left = mid + 1;
    }

    if (left < vOligonucleotideList.size() && vOligonucleotideList[left].dMass + (Param::g_vVarModification[vOligonucleotideList[left].iVarModificationIndex]).varModMassRange.dStart >= dMass) 
    {
        return left;
    }

    return 0;
}

int CalculateUtils::BinarySearchOligonucleotideByHighestMass(const vector<Oligonucleotide>& vOligonucleotideList, double dMass) //should sort by highest mass first
{
    int left = 0;
    int right = vOligonucleotideList.size() - 1;

    if (right < 0) {
        return -1;
    }

    while (left <= right) {
        int mid = left + (right - left) / 2;
        if (vOligonucleotideList[mid].dMass + (Param::g_vVarModification[vOligonucleotideList[mid].iVarModificationIndex]).varModMassRange.dEnd > dMass)
            right = mid - 1;
        else
            left = mid + 1;
    }

    if (right < vOligonucleotideList.size() && right >= 0 && vOligonucleotideList[right].dMass + (Param::g_vVarModification[vOligonucleotideList[right].iVarModificationIndex]).varModMassRange.dEnd < dMass) {
        return right;
    }

    return 0;
}


//////////////////////////////////////////////////////////////////////
// Implementations for CombinatoricsUtils
//////////////////////////////////////////////////////////////////////

int** BINOM_COEF; // array[N][K], partial Pascal's Triangle
int      N = -1;
int      K = -1;

void CombinatoricsUtils::initBinomialCoefficients(const int n, const int k)
{
    N = n;
    K = k;
    BINOM_COEF = new int* [n + 1];
    for (int i = 0; i <= n; ++i)
    {
        const auto coeffs = new int[k + 1];
        if (i == 0)
        {
            coeffs[0] = 1;
        }
        else
        {

            for (int j = 0; j <= i && j <= k; ++j)
            {
                if (j == 0 || j == i)
                {
                    coeffs[j] = 1;
                    continue;
                }
                int* prevRow = BINOM_COEF[i - 1];
                coeffs[j] = prevRow[j - 1] + prevRow[j];
            }
        }
        BINOM_COEF[i] = coeffs;
    }
}

// Copied from org.apache.commons.math3.util.CombinatoricsUtils.java
// Returns a two dimensional array [count][r]
int** CombinatoricsUtils::makeCombinations(int n, int r, int count)
{
    if (r == 0)
    {
        return new int* [0];
    }

    if (r == n)
    {
        int* combination = new int[n];
        for (int i = 0; i < n; ++i)
        {
            combination[i] = i;
        }
        int** combinations = new int* [1];
        combinations[0] = combination;
        return combinations;
    }

    int k = r;
    int j = k; // Set up invariant: j is smallest index such that c[j + 1] > j
    int* c = new int[k + 3];

    c[0] = 0;
    for (int i = 1; i <= k; ++i)
    {
        c[i] = i - 1;
    }
    // Initialize sentinels
    c[k + 1] = n;
    c[k + 2] = 0;

    bool more = true;

    int** combinations = new int* [count];
    int idx = 0;
    while (more)
    {
        int* ret = new int[k];
        // Copy return value (prepared by last activation)
        // System.arraycopy(c, 1, ret, 0, k);
        for (int i = 0; i < k; ++i)
        {
            ret[i] = c[i + 1];
        }

        // Prepare next iteration
        // T2 and T6 loop
        int x = 0;
        if (j > 0)
        {
            x = j;
            c[j] = x;
            j--;
            combinations[idx] = ret;
            idx++;
            continue;
            // return ret;
        }
        // T3
        if (c[1] + 1 < c[2])
        {
            c[1]++;
            combinations[idx] = ret;
            idx++;
            continue;
            // return ret;
        }
        else
        {
            j = 2;
        }
        // T4
        bool stepDone = false;
        while (!stepDone)
        {
            c[j - 1] = j - 2;
            x = c[j] + 1;
            if (x == c[j + 1])
            {
                j++;
            }
            else
            {
                stepDone = true;
            }
        }
        // T5
        if (j > k)
        {
            more = false;
            combinations[idx] = ret;
            idx++;
            continue;
            // return ret;
        }
        // T6
        c[j] = x;
        j--;
        combinations[idx] = ret;
        idx++;
        // return ret;
    }
    delete[] c;
    return combinations;
}

int CombinatoricsUtils::nChooseK(const int n, const int k)
{
    if (n == k || k == 0)
        return 1;

    if (k == 1 || k == n - 1)
        return n;

    if (k > n / 2)
        return nChooseK(n, n - k);

    if (n <= N && k <= K)
    {
        return BINOM_COEF[n][k];
    }

    // https://en.wikipedia.org/wiki/Binomial_coefficient#Identities_involving_binomial_coefficients
    // (n, k) = n / k * (n - 1, k - 1)
    int answer = n - k + 1; // base (n - k + 1, 1)
    int previous = answer;
    for (int i = 1; i < k; ++i)
    {
        answer = previous * (n - k + 1 + i) / (i + 1);
        previous = answer;
    }

    return answer;
}

int CombinatoricsUtils::getCombinationCount(int n, int k)
{
    int total = 0;
    if (k > n)
        k = n;
    for (; k >= 1; k--)
    {
        total += nChooseK(n, k);
    }
    return total;
}









