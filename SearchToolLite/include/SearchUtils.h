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

#include <vector>
#include <string>
#include <map>

class SearchUtils
{
public:
	SearchUtils();
   ~SearchUtils();

   bool PrecursorFilter(ThreadPool* _tp, Query& query, const vector<Oligonucleotide>& vOligo);

   bool XcorrScoring(ThreadPool* _tp, Query& query);
   double XcorrScoring(ThreadPool* _tp, const ExpSpectrum& expSpectrum, const vector<OligoFragment>& vFragments);

   bool EntropyScoring(ThreadPool* _tp, Query& query);
   double EntropyScoring(ThreadPool* _tp, const ExpSpectrum& expSpectrum, const vector<OligoFragment>& vFragments);

   bool IonMatchScoring(ThreadPool* _tp, Query& query);
   double IonMatchScoring(ThreadPool* _tp, const ExpSpectrum& expSpectrum, const vector<OligoFragment>& vFragments);

   static SearchUtils& getInstance()
   {
      static SearchUtils instance;
      return instance;
   }

   double dTimeFragment;
   double dTimeXcorrCalculate;

private:
   // Core search functions
   double* CalTheoreticalSpectrum(const vector<OligoFragment>& vFragments, int iArraySize);

   double* CalExperimentalSpectrum(const ExpSpectrum& expSpectrum);

   // {peak_index, fragment_index} 
   vector<pair<int, int>> FragmentMatchCount(const ExpSpectrum& expSpectrum, const vector<OligoFragment>& vFragmentList);

   std::map<OligonucleotideID, vector<OligoFragment>> mOligoIDFragments;
};

// to speed up the precursor mass search, we use a index to store the precursor mass
class PrecursorMassIndexUtils
{
public:
    struct PrecursorMassIndex
    {
        double dMass;
        OligonucleotideID oligoID;

        PrecursorMassIndex(double dMass, int iOligoIndex, int iPositionIndex, int iVarModFormulaIndex, int iVarModIndex)
        {
            this->dMass = dMass;
            this->oligoID = OligonucleotideID(iOligoIndex, iPositionIndex, iVarModFormulaIndex, iVarModIndex);
        }
        
        PrecursorMassIndex(const PrecursorMassIndex& other)
        {
            this->dMass = other.dMass;
            this->oligoID = other.oligoID;
        }

        PrecursorMassIndex& operator=(const PrecursorMassIndex& other)
        {
            if (this != &other)
            {
                dMass = other.dMass;
                oligoID = other.oligoID;
            }
            return *this;
        }
    };


    PrecursorMassIndexUtils();
    ~PrecursorMassIndexUtils();

    bool GeneratePrecursorMassIndex(const vector<Oligonucleotide>& vOligo);

    vector<vector<PrecursorMassIndex>> GetIndices(const DoubleRange& massRange);

    static PrecursorMassIndexUtils& getInstance()
    {
        static PrecursorMassIndexUtils instance;
        return instance;
    }

private:
    vector<vector<PrecursorMassIndex>> vvPrecursorMassIndices;

};

// to speed up the fragment m/z match search, we use a index to store the fragment m/z
class FragmentMZIndexUtils
{
    struct MZIndex
    {
        double dMZ;
        const OligoFragment& fragment;

        MZIndex(double dMZ, const OligoFragment& fragment) : fragment(fragment)
        {
            this->dMZ = dMZ;
        }
    };

public:
    FragmentMZIndexUtils();
    ~FragmentMZIndexUtils();

    bool GenerateFragmentMZIndex(const vector<OligoFragment>& vFragments);

    std::vector<std::vector<MZIndex>> GetIndices(const DoubleRange& mzRange);

    static FragmentMZIndexUtils& getInstance()
    {
        static FragmentMZIndexUtils instance;
        return instance;
    }

private:
    std::vector<std::vector<MZIndex>> vvFragmentMZIndices;
};

#endif
