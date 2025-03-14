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


#include "SearchUtils.h"
#include "MassUtils.h"
#include "CalculateUtils.h"
#include "OligoProcessUtils.h"
#include <vector>


SearchUtils::SearchUtils()
{
}

SearchUtils::~SearchUtils()
{  
}

bool SearchUtils::PrecursorFilter(ThreadPool* _tp, Query& query, const vector<Oligonucleotide>& vOligo)
{
   for (const auto& range : query.expSpectrum.vPreMassMatchRange)
   {
      vector<vector<PrecursorMassIndexUtils::PrecursorMassIndex>> vvPrecursorMassIndices = PrecursorMassIndexUtils::getInstance().GetIndices(range);
      for (const auto& vPrecursorMassIndex : vvPrecursorMassIndices)
      {
         for (const auto& precursorMassIndex : vPrecursorMassIndex)
         {
            const Oligonucleotide& oligonucleotide = vOligo[precursorMassIndex.oligoID.iOligoIndex];
            const VarModFormula& varModFormula = Param::getInstance().vVarModFormulas[precursorMassIndex.oligoID.iVarModFormulaIndex];
            Query::OligoEntry oligoEntry(precursorMassIndex.oligoID.iOligoIndex, precursorMassIndex.oligoID.iPositionIndex, precursorMassIndex.oligoID.iVarModFormulaIndex, precursorMassIndex.oligoID.iVarModIndex);
            query.vOligoEntries.push_back(oligoEntry);
         }
      }
   }

   return true;
}

bool SearchUtils::XcorrScoring(ThreadPool* _tp, Query& query)
{
   if (query.vOligoEntries.size() == 0)
       return true;

   for (auto& oligoEntry: query.vOligoEntries)
   {
      auto startTime = std::chrono::high_resolution_clock::now();
      vector<OligoFragment> vFragments = OligoProcessUtils::getInstance().GenerateFragments(Param::getInstance().vOligos[oligoEntry.oligoID.iOligoIndex], query.expSpectrum.iPreCharge, oligoEntry.oligoID);
      auto endTime = std::chrono::high_resolution_clock::now();
      this->dTimeFragment += std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count() / 1000.0;

      startTime = std::chrono::high_resolution_clock::now();
      oligoEntry.scores.xCorr = XcorrScoring(_tp, query.expSpectrum, vFragments);
      endTime = std::chrono::high_resolution_clock::now();
      this->dTimeXcorrCalculate += std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count() / 1000.0;
   }

   std::sort(query.vOligoEntries.begin(), query.vOligoEntries.end(), [](const Query::OligoEntry& a, const Query::OligoEntry& b) {
       return a.scores.xCorr > b.scores.xCorr; // Sort in descending order
   });

   return true;
}

double SearchUtils::XcorrScoring(ThreadPool* _tp, const ExpSpectrum& expSpectrum, const vector<OligoFragment>& vFragments)
{
   if (vFragments.size() == 0)
      return true;
   
   double* theoreticalSpectrum = CalTheoreticalSpectrum(vFragments, expSpectrum.iArraySize);
   double* experimentalSpectrum = CalExperimentalSpectrum(expSpectrum);
   double dXcorr = CalculateUtils::getInstance().FastXcorrCalculate(theoreticalSpectrum, experimentalSpectrum, expSpectrum.iArraySize, Param::getInstance().options.iXcorrProcessingOffset);
   if (dXcorr < CommonValues::dFloatZero)
      dXcorr = 0;
   delete[] theoreticalSpectrum;
   delete[] experimentalSpectrum;
   return dXcorr;
}

bool SearchUtils::EntropyScoring(ThreadPool* _tp, Query& query)
{
   if (query.vOligoEntries.size() == 0)
      return true;

   // generate theoretical peaks
   vector<ExpSpectrum::Peak> vTheoreticalPeaks;
   for (auto& oligoEntry : query.vOligoEntries)
   {
      vector<OligoFragment> vFragments = OligoProcessUtils::getInstance().GenerateFragments(Param::getInstance().vOligos[oligoEntry.oligoID.iOligoIndex], query.expSpectrum.iPreCharge, oligoEntry.oligoID);
      for (const auto& fragment : vFragments)
      {
         if (fragment.iCharge != 0)
         {
            double dMZ = MassUtils::getInstance().MassToMZ(fragment.dMass, fragment.iCharge, Param::getInstance().options.sAdductIon);
            vTheoreticalPeaks.push_back(ExpSpectrum::Peak(dMZ, 1.0 / vFragments.size()));
         }
         else
         {
            int iCharge;
            if (Param::getInstance().options.bPositiveCharge)
               iCharge = 1;
            else
               iCharge = -1;
            
            double dMZ = MassUtils::getInstance().MassToMZ(fragment.dMass, iCharge, Param::getInstance().options.sAdductIon);
            vTheoreticalPeaks.push_back(ExpSpectrum::Peak(dMZ, 1.0 / vFragments.size()));
         }
      }
      // get experimental peaks
      vector<ExpSpectrum::Peak> vExperimentalPeaks = query.expSpectrum.spectrum;
      double dHighestIntensity = vExperimentalPeaks[query.expSpectrum.iHighestIonIndex].dIntensity;
      for (auto& peak : vExperimentalPeaks)
         peak.dIntensity /= dHighestIntensity;  
      
      // merge peaks
      vector<ExpSpectrum::Peak> vMergedPeaks;
      {
         DoubleRange tolerance = Param::getInstance().options.fragmentTolerance;
         int iIndex1 = 0;
         int iIndex2 = 0;
         while (iIndex1 < vTheoreticalPeaks.size() && iIndex2 < vExperimentalPeaks.size())
         {
            if (vMergedPeaks.size() > 0)
            {
               double dTheoreticalPPM = (vTheoreticalPeaks[iIndex1].dMZ - vMergedPeaks.back().dMZ) / vMergedPeaks.back().dMZ;
               double dExperimentalPPM = (vExperimentalPeaks[iIndex2].dMZ - vMergedPeaks.back().dMZ) / vMergedPeaks.back().dMZ;
               if (tolerance.inRange(dTheoreticalPPM))
               {
                  vMergedPeaks.back().dIntensity += vTheoreticalPeaks[iIndex1].dIntensity;
                  iIndex1++;
               }
               else if (tolerance.inRange(dExperimentalPPM))
               {
                  vMergedPeaks.back().dIntensity += vExperimentalPeaks[iIndex2].dIntensity;
                  iIndex2++;
               }
               else
               {
                  if (vTheoreticalPeaks[iIndex1].dMZ < vExperimentalPeaks[iIndex2].dMZ)
                     vMergedPeaks.push_back(vTheoreticalPeaks[iIndex1++]);
                  else
                     vMergedPeaks.push_back(vExperimentalPeaks[iIndex2++]);
               }
            }
            else
               if (vTheoreticalPeaks[iIndex1].dMZ < vExperimentalPeaks[iIndex2].dMZ)
                  vMergedPeaks.push_back(vTheoreticalPeaks[iIndex1++]);
               else
                  vMergedPeaks.push_back(vExperimentalPeaks[iIndex2++]);
         }
      }

      // calculate entropy
      auto calculateEntropy = [](const vector<ExpSpectrum::Peak>& vPeaks)
      {
         double dEntropy = 0;
         for (const auto& peak : vPeaks)
            dEntropy += -1 * peak.dIntensity * log(peak.dIntensity);
         return dEntropy;
      };
      double dTheoreticalEntropy = calculateEntropy(vTheoreticalPeaks);
      double dExperimentalEntropy = calculateEntropy(vExperimentalPeaks);
      double dMergedEntropy = calculateEntropy(vMergedPeaks);
      oligoEntry.scores.entropy = 1 - (2 * dMergedEntropy - dTheoreticalEntropy - dExperimentalEntropy) / log(4);
      
   }

   std::sort(query.vOligoEntries.begin(), query.vOligoEntries.end(), [](const Query::OligoEntry& a, const Query::OligoEntry& b) {
       return a.scores.entropy > b.scores.entropy; // Sort in descending order
   });

   return true;
}

double SearchUtils::EntropyScoring(ThreadPool* _tp, const ExpSpectrum& expSpectrum, const vector<OligoFragment>& vFragments)
{
   return 0;
}

bool SearchUtils::IonMatchScoring(ThreadPool* _tp, Query& query)
{
   if (query.vOligoEntries.size() == 0)
       return true;

   for (auto& oligoEntry: query.vOligoEntries)
   {
      vector<OligoFragment> vFragments = OligoProcessUtils::getInstance().GenerateFragments(Param::getInstance().vOligos[oligoEntry.oligoID.iOligoIndex], query.expSpectrum.iPreCharge, oligoEntry.oligoID);
      vector<pair<int, int>> vMatchInfo = SearchUtils::getInstance().FragmentMatchCount(query.expSpectrum, vFragments);
      oligoEntry.scores.iMatchIonCount = vMatchInfo.size();
      oligoEntry.scores.iTotalIonCount = vFragments.size();
      oligoEntry.scores.dSp = SearchUtils::getInstance().IonMatchScoring(_tp, query.expSpectrum, vFragments);
   }

   std::sort(query.vOligoEntries.begin(), query.vOligoEntries.end(), [](const Query::OligoEntry& a, const Query::OligoEntry& b) {
       return a.scores.dSp > b.scores.dSp; // Sort in descending order
   });

   return true;
}

double SearchUtils::IonMatchScoring(ThreadPool* _tp, const ExpSpectrum& expSpectrum, const vector<OligoFragment>& vFragments)
{
   if (vFragments.size() == 0)
      return true;

   vector<pair<int, int>> vMatchInfo = SearchUtils::getInstance().FragmentMatchCount(expSpectrum, vFragments);

   // consecutive ion match
   map<IonSpecies::IonID, vector<int>> mIonSpeciesFragmentIndices;
   for (const auto& matchInfo : vMatchInfo)
   {
      mIonSpeciesFragmentIndices[vFragments[matchInfo.second].ionType.ionType].push_back(matchInfo.second);
   }

   // calculate Sp
   double dSp = 0;
   double dMatchIntensity = 0;
   int iTotalIons = vFragments.size();
   int iMatchIons = 0;
   for (auto& pair : vMatchInfo)
   {
      double dIntensity = expSpectrum.spectrum[pair.first].dIntensity;
      dMatchIntensity += dIntensity;
      iMatchIons++;
   }
   dSp = dMatchIntensity / expSpectrum.dTotalIntensity * ((double)(iMatchIons)/(double)(iTotalIons));
   return dSp;
}

double* SearchUtils::CalTheoreticalSpectrum(const vector<OligoFragment>& vFragments, int iArraySize)
{
   if (iArraySize <= 0)
      return nullptr;

   double* spectrum = new double[iArraySize];
   for (int i = 0; i < iArraySize; i++)
      spectrum[i] = 0;

   for (const auto& fragment: vFragments)
   {
      if (fragment.iCharge != 0)
      {
         double dMZ = MassUtils::getInstance().MassToMZ(fragment.dMass, fragment.iCharge, Param::getInstance().options.sAdductIon);
         int iBinIndex = MassUtils::getInstance().CalBinIndex(dMZ);
         if (iBinIndex <= iArraySize - 1)
            spectrum[iBinIndex] += 1.0;
      }
      else
      {
         int iCharge;
         if (Param::getInstance().options.bPositiveCharge)
            iCharge = 1;
         else
            iCharge = -1;
         
         double dMZ = MassUtils::getInstance().MassToMZ(fragment.dMass, iCharge, Param::getInstance().options.sAdductIon);
         int iBinIndex = MassUtils::getInstance().CalBinIndex(dMZ);
         if (iBinIndex <= iArraySize - 1)
            spectrum[iBinIndex] += 1.0;
      }
   }

   return spectrum;
}

double* SearchUtils::CalExperimentalSpectrum(const ExpSpectrum& expSpectrum)
{
   if (expSpectrum.iArraySize <= 0)
      return nullptr;

   double* spectrum = new double[expSpectrum.iArraySize];
   for (int i = 0; i < expSpectrum.iArraySize; i++)
      spectrum[i] = 0;

   const vector<ExpSpectrum::Peak>& vPeaks = expSpectrum.spectrum;
   double dMaxIntensity = vPeaks[expSpectrum.iHighestIonIndex].dIntensity;
   for (const auto& peak : vPeaks)
   {
      double dMZ = peak.dMZ;
      int iBinIndex = MassUtils::getInstance().CalBinIndex(dMZ);
      if (iBinIndex < expSpectrum.iArraySize)
          spectrum[iBinIndex] += peak.dIntensity / dMaxIntensity;
   }
   
   return spectrum;
}

vector<pair<int, int>> SearchUtils::FragmentMatchCount(const ExpSpectrum& expSpectrum, const vector<OligoFragment>& vFragmentList)
{
   vector<pair<int, int>> vMatchInfo;
   const vector<ExpSpectrum::Peak>& vPeaks = expSpectrum.spectrum;
   for (size_t iPeakIndex = 0; iPeakIndex < vPeaks.size(); iPeakIndex++)
   {
      const ExpSpectrum::Peak& peak = vPeaks[iPeakIndex];
      if (peak.dIntensity < CommonValues::dFloatZero)
         continue;

      for (size_t iFragmentIndex = 0; iFragmentIndex < vFragmentList.size(); iFragmentIndex++)
      {
         const OligoFragment& fragment = vFragmentList[iFragmentIndex];
         double dMass = fragment.dMass;
         int iCharge = fragment.iCharge;
         double dMZ = MassUtils::getInstance().MassToMZ(dMass, iCharge, Param::getInstance().options.sAdductIon);

         DoubleRange tolerance = Param::getInstance().options.fragmentTolerance;
         tolerance *= ((double)dMZ / (double)1e6);
         if (tolerance.inRange(peak.dMZ - dMZ))
         {
            vMatchInfo.push_back({iPeakIndex, iFragmentIndex});
            break;   
         }
         
      }
   }
   return vMatchInfo;
}

// ==============================
// PrecursorMassIndexUtils
// ==============================

PrecursorMassIndexUtils::PrecursorMassIndexUtils()
{
   int iBinNumber = int(Param::getInstance().options.digestMassRange.dEnd) + 1;
   vvPrecursorMassIndices = vector<vector<PrecursorMassIndex>>(iBinNumber);
}

PrecursorMassIndexUtils::~PrecursorMassIndexUtils()
{
}

bool PrecursorMassIndexUtils::GeneratePrecursorMassIndex(const vector<Oligonucleotide>& vOligo)
{
   for (int iOligoIndex = 0; iOligoIndex < vOligo.size(); iOligoIndex++)
   {
      const Oligonucleotide& oligonucleotide = vOligo[iOligoIndex];

      for (const auto& iVarModFormulaIndex : oligonucleotide.vVarModFormulaIndices)
      {
         const VarModFormula& varModFormula = Param::getInstance().vVarModFormulas[iVarModFormulaIndex];
         for (int iVarModIndex = 0; iVarModIndex < varModFormula.vVarModMassList.size(); iVarModIndex++)
         {
            double dVarModMass = varModFormula.vVarModMassList[iVarModIndex];
            OligonucleotideID oligoID(iOligoIndex, 0, iVarModFormulaIndex, iVarModIndex);
            double dOligoMass = MassUtils::getInstance().CalOligoMass(oligonucleotide, oligoID, true, true, true) + dVarModMass;
            int iBinIndex = int(dOligoMass);
            if (iBinIndex < vvPrecursorMassIndices.size())
               for (int iPositionIndex = 0; iPositionIndex < oligonucleotide.mVarModPosition.at(iVarModFormulaIndex).size(); iPositionIndex++)
                  vvPrecursorMassIndices[iBinIndex].push_back(PrecursorMassIndex(dOligoMass, iOligoIndex, iPositionIndex, iVarModFormulaIndex, iVarModIndex));
         }
      }
   }
   return true;
}

vector<vector<PrecursorMassIndexUtils::PrecursorMassIndex>> PrecursorMassIndexUtils::GetIndices(const DoubleRange& massRange)
{
   vector<vector<PrecursorMassIndex>> vvResult;
   int iStartIndex = int(massRange.dStart);
   int iEndIndex = int(massRange.dEnd);
   if (iStartIndex < 0)
      iStartIndex = 0;
   if (iEndIndex >= vvPrecursorMassIndices.size())
      iEndIndex = vvPrecursorMassIndices.size() - 1;

   for (size_t iIndex = iStartIndex; iIndex <= iEndIndex; iIndex++)
   {
      vector<PrecursorMassIndex> vResult;
      for (PrecursorMassIndex& precursorMassIndex : vvPrecursorMassIndices[iIndex])
      {
         if (massRange.inRange(precursorMassIndex.dMass))
            vResult.push_back(precursorMassIndex);
      }
      vvResult.push_back(vResult);
   }
   return vvResult;
}

// ==============================
// FragmentMZIndexUtils
// ==============================

FragmentMZIndexUtils::FragmentMZIndexUtils()
{
   int iBinNumber = int(Param::getInstance().options.fragmentMassRange.dEnd) + 1;
   vvFragmentMZIndices = vector<vector<MZIndex>>(iBinNumber);
}

FragmentMZIndexUtils::~FragmentMZIndexUtils()
{
}  

bool FragmentMZIndexUtils::GenerateFragmentMZIndex(const vector<OligoFragment>& vFragments)
{
   for (const auto& fragment : vFragments)
   {
      if (fragment.iCharge != 0)
      {
         double dMZ = MassUtils::getInstance().MassToMZ(fragment.dMass, fragment.iCharge, Param::getInstance().options.sAdductIon);
         int iBinIndex = MassUtils::getInstance().CalBinIndex(dMZ);
         if (iBinIndex < vvFragmentMZIndices.size())
            vvFragmentMZIndices[iBinIndex].push_back(MZIndex(dMZ, fragment));
      }
      else
      {
         int iCharge;
         if (Param::getInstance().options.bPositiveCharge)
            iCharge = 1;
         else
            iCharge = -1;
         double dMZ = MassUtils::getInstance().MassToMZ(fragment.dMass, iCharge, Param::getInstance().options.sAdductIon);
         int iBinIndex = MassUtils::getInstance().CalBinIndex(dMZ);
         if (iBinIndex < vvFragmentMZIndices.size())
            vvFragmentMZIndices[iBinIndex].push_back(MZIndex(dMZ, fragment));
      }
   }
   return true;
}

std::vector<std::vector<FragmentMZIndexUtils::MZIndex>> FragmentMZIndexUtils::GetIndices(const DoubleRange& dMZRange)
{
   return std::vector<std::vector<FragmentMZIndexUtils::MZIndex>>();
}

