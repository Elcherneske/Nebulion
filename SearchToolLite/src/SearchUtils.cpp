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

bool SearchUtils::RoughScoring(ThreadPool* _tp, Query& query)
{
   if (query.vOligoEntries.size() == 0)
       return true;

   for (auto& oligoEntry: query.vOligoEntries)
   {
       vector<OligoFragment> vFragments;
      if (mOligoIDFragments.find(oligoEntry.oligoID) != mOligoIDFragments.end())
          vFragments = mOligoIDFragments[oligoEntry.oligoID];
      else
      {
          vFragments = OligoProcessUtils::getInstance().GenerateFragments(Param::getInstance().vOligos[oligoEntry.oligoID.iOligoIndex], query.expSpectrum.iPreCharge, oligoEntry.oligoID);
          mOligoIDFragments[oligoEntry.oligoID] = vFragments;
      }

      int* theoreticalSpectrum = CalTheoreticalSpectrum(vFragments, query.expSpectrum.iArraySize);
      int* experimentalSpectrum = CalExperimentalSpectrum(query.expSpectrum);
      double dXcorr = CalculateUtils::getInstance().FastXcorrCalculate(theoreticalSpectrum, experimentalSpectrum, query.expSpectrum.iArraySize, Param::getInstance().options.iXcorrProcessingOffset);
      oligoEntry.scores.xCorr = dXcorr;
      if (oligoEntry.scores.xCorr < CommonValues::dFloatZero)
         oligoEntry.scores.xCorr = 0;
      delete[] theoreticalSpectrum;
      delete[] experimentalSpectrum;
   }

   std::sort(query.vOligoEntries.begin(), query.vOligoEntries.end(), [](const Query::OligoEntry& a, const Query::OligoEntry& b) {
       return a.scores.xCorr > b.scores.xCorr; // Sort in descending order
   });
   return true;
}

// bool SearchUtils::EntryScoring(const ExpSpectrum& expSpectrum, const Oligonucleotide& oligonucleotide, int iWhichVarModCombination, Scores& result)
// {
//    bool bSuccessed = true;
//    if (abs(expSpectrum.iPreCharge) == 0)
//    {
//        cout << "Error: in scoring procedure, the spectrum precursor state is 0" << endl;
//        return false;
//    }
//    int iCharge = expSpectrum.iPreCharge;

//    //generate fragments
//    vector<OligonucleotideFragment> vfragmentList;
//    for (IonSeries::Ions ionType : Param::g_staticParams.ionInformation.vSelectedIonSeries)
//    {
//       vector<pair<int, Modification>> varModCombination = Param::g_vVarModification[oligonucleotide.iVarModificationIndex].varModificationList[iWhichVarModCombination];
//       vector<OligonucleotideFragment> vTmpFragmentList = OligoProcessUtils::GenerateOligonucleotideFragment(oligonucleotide, ionType, iCharge, varModCombination);
//       vfragmentList.insert(vfragmentList.end(), vTmpFragmentList.begin(), vTmpFragmentList.end());
//    }
//    result.totalIons = vfragmentList.size();

//    //generate theoretical spectrum and expeimental spectrum
//    vector<double> theoreticalSpectrum;
//    vector<double> experimentalSpectrum;

//    bSuccessed = GenerateTheoreticalSpectrum(vfragmentList, expSpectrum.iArraySize, iCharge, theoreticalSpectrum);
//    if (!bSuccessed)
//        return false;
//    bSuccessed = GenerateExperimentalSpectrum(expSpectrum, experimentalSpectrum);
//    if (!bSuccessed)
//        return false;

//    // calculate Xcorr
//    double dXcorr;
//    bSuccessed = CalculateUtils::FastXcorrCalc(theoreticalSpectrum, experimentalSpectrum, Param::g_staticParams.iXcorrProcessingOffset, dXcorr);
//    if (!bSuccessed)
//        return false;
//    result.xCorr = dXcorr;
//    if (result.xCorr < CommonValues::dFloatZero)
//       result.xCorr = 0;

//    //cout match fragment
//    vector<pair<int, int>> vIonSpectrumMatch = SearchUtils::FragmentMatchCount(expSpectrum, vfragmentList); //{peak, fragment}
//    result.matchedIons = vIonSpectrumMatch.size();

//    //calculate Sp bonus
//    int iConsec = 0;
//    for (int i = 0; i < vIonSpectrumMatch.size(); i++)
//    {
//       for (int j = i + 1; j < vIonSpectrumMatch.size(); j++)
//       {
//          if (vfragmentList[vIonSpectrumMatch[i].second].ionType == vfragmentList[vIonSpectrumMatch[j].second].ionType
//          && vfragmentList[vIonSpectrumMatch[i].second].sSequence.size() == vfragmentList[vIonSpectrumMatch[j].second].sSequence.size() - 1)
//             iConsec++;
//       }
//    }

//    //calculate Sp
//    double dMatchIntensity = 0.0;
//    Spectrum spec = Spectrum(expSpectrum.spectrum);
//    double dHighestIntensity = spec.at(expSpectrum.iHighestIonIndex).intensity;
//    double dTotalIntensity = 0.0;
//    for (auto& pair : vIonSpectrumMatch)
//    {
//       double dIntensity = spec.at(pair.first).intensity;
//       dMatchIntensity += dIntensity;
//    }

//    result.dSp = dMatchIntensity/ expSpectrum.dTotalIntensity * ((double)(1+iConsec)/(double)(result.totalIons));
//    if (result.dSp < CommonValues::dFloatZero)
//       result.dSp = 0;
//    return true;
// }

int* SearchUtils::CalTheoreticalSpectrum(const vector<OligoFragment>& vFragments, int iArraySize)
{
   if (iArraySize <= 0)
      return nullptr;

   int* spectrum = new int[iArraySize];
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

int* SearchUtils::CalExperimentalSpectrum(const ExpSpectrum& expSpectrum)
{
   if (expSpectrum.iArraySize <= 0)
      return nullptr;

   int* spectrum = new int[expSpectrum.iArraySize];
   for (int i = 0; i < expSpectrum.iArraySize; i++)
      spectrum[i] = 0;

   const vector<ExpSpectrum::Peak>& vPeaks = expSpectrum.spectrum;
   double dMaxIntensity = vPeaks[expSpectrum.iHighestIonIndex].dIntensity;
   for (const auto& peak : vPeaks)
   {
      double dMZ = peak.dMZ;
      int iBinIndex = MassUtils::getInstance().CalBinIndex(dMZ);
      if (iBinIndex < expSpectrum.iArraySize)
          spectrum[iBinIndex] += (double)peak.dIntensity / dMaxIntensity;
   }
   
   return spectrum;
}

// vector<pair<int, int>> SearchUtils::FragmentMatchCount(const ExpSpectrum& expSpectrum, const vector<OligonucleotideFragment>& vFragmentList)
// {
//    int iMatchCount = 0;
//    vector<pair<int, int>> vMatchIndex;
//    Spectrum spectrum = Spectrum(expSpectrum.spectrum);
//    for (int i = 0; i < spectrum.size(); i++)
//    {
//       Peak_T peak = spectrum.at(i);
//       if (peak.intensity < CommonValues::dFloatZero)
//          continue;
//       double dExpMZ = peak.mz;

//       for (int ii = 0; ii < vFragmentList.size(); ii++)
//       {
//          const OligonucleotideFragment& fragment = vFragmentList[ii];
//          double dMass = fragment.dMass;
//          int iCharge = fragment.iCharge;

//          if (Param::g_staticParams.tolerances.iFragmentToleranceUnits == 0) // amu
//          {
//             if (Param::g_staticParams.tolerances.iFragmentToleranceType == 1)  // precursor m/z tolerance
//             {
//                double dTheMZ = MassUtils::MassToMZ(dMass, iCharge);
//                if (DoubleRange::inRange(Param::g_staticParams.tolerances.fragmentTolerance, dExpMZ - dTheMZ))
//                {
//                   iMatchCount ++;
//                   vMatchIndex.push_back({i, ii});
//                   break;
//                }
//             }
//             else
//             {
//                double dExpMass = MassUtils::MZToMass(dExpMZ, iCharge);
//                if (DoubleRange::inRange(Param::g_staticParams.tolerances.fragmentTolerance, dExpMass - dMass))
//                {
//                   iMatchCount ++;
//                   vMatchIndex.push_back({i, ii});
//                   break;
//                }
//             }
//          }
//          else if (Param::g_staticParams.tolerances.iFragmentToleranceUnits == 1) // mmu
//          {
//             DoubleRange tolerance = Param::g_staticParams.tolerances.fragmentTolerance;
//             tolerance *= 0.001;
//             if (Param::g_staticParams.tolerances.iFragmentToleranceType == 1)  // precursor m/z tolerance
//             {
//                double dTheMZ = MassUtils::MassToMZ(dMass, iCharge);
//                if (DoubleRange::inRange(tolerance, dExpMZ - dTheMZ))
//                {
//                   iMatchCount ++;
//                   vMatchIndex.push_back({i, ii});
//                   break;
//                }
//             }
//             else
//             {
//                double dExpMass = MassUtils::MZToMass(dExpMZ, iCharge);
//                if (DoubleRange::inRange(tolerance, dExpMass - dMass))
//                {
//                   iMatchCount ++;
//                   vMatchIndex.push_back({i, ii});
//                   break;
//                }
//             }
//          }
//          else // ppm
//          {
//             double dTheMZ = MassUtils::MassToMZ(dMass, iCharge);

//             DoubleRange tolerance = Param::g_staticParams.tolerances.fragmentTolerance;
//             tolerance *= ((double)dTheMZ / (double)1e6);
//             if (DoubleRange::inRange(tolerance, dExpMZ - dTheMZ))
//             {
//                iMatchCount ++;
//                vMatchIndex.push_back({i, ii});
//                break;
//             }
//          }
//       }
//    }
//    return vMatchIndex;
// }

// PrecursorMassIndexUtils

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

// FragmentMZIndexUtils
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

