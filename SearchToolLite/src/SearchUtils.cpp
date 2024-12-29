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


SearchUtils::SearchUtils()
{
}

SearchUtils::~SearchUtils()
{  
}

bool SearchUtils::MS1Search(Query& query, vector<Oligonucleotide>& vOligonucleotideList)
{
    if (query.expSpectrum == nullptr)
        return false;

    int index = CalculateUtils::BinarySearchOligonucleotideByHighestMass(vOligonucleotideList, query.expSpectrum->preMassSearchRange.dStart);
    if (index >= 0 && index < vOligonucleotideList.size())
    {
        for (int i = index; i < vOligonucleotideList.size(); i++)
        {
            Oligonucleotide& oligonucleotide = vOligonucleotideList[i];
            const VarModification& varMod = Param::g_vVarModification[oligonucleotide.iVarModificationIndex];

            if (oligonucleotide.dMass + varMod.varModMassRange.dStart > query.expSpectrum->preMassSearchRange.dEnd)
                break;

            for (int i = 0; i < varMod.vVarModificationMassList.size(); i++)
            {
                double dVarModMass = varMod.vVarModificationMassList[i];
                for (DoubleRange& searchRange : query.expSpectrum->vPreMassSearchRange)
                {
                    if (searchRange.inRange(oligonucleotide.dMass + dVarModMass))
                    {
                        query.pushBackResult(oligonucleotide, i);
                        break;
                    }
                }
            }
        }
    }
    else
        return false;

    return true;
}

bool SearchUtils::QueryScoring(ThreadPool* _tp, Query& query, const vector<Oligonucleotide>& vOligonucleotideList)
{
   if(query.expSpectrum == nullptr)
      return false;

   if (query.vResults.size() == 0)
       return true;

   // do MS2 match/score
   for (auto& resultItem: query.vResults)
   {
      SearchUtils::EntryScoring(*(query.expSpectrum), resultItem.oligo, resultItem.iWhichVarModCombination, resultItem.scores);
   }

   return true;
}

bool SearchUtils::EntryScoring(const ExpSpectrum& expSpectrum, const Oligonucleotide& oligonucleotide, int iWhichVarModCombination, Scores& result)
{
   bool bSuccessed = true;
   if (abs(expSpectrum.iPreCharge) == 0)
   {
       cout << "Error: in scoring procedure, the spectrum precursor state is 0" << endl;
       return false;
   }
   int iCharge = expSpectrum.iPreCharge;

   //generate fragments
   vector<OligonucleotideFragment> vfragmentList;
   for (IonSeries::Ions ionType : Param::g_staticParams.ionInformation.vSelectedIonSeries)
   {
      vector<pair<int, Modification>> varModCombination = Param::g_vVarModification[oligonucleotide.iVarModificationIndex].varModificationList[iWhichVarModCombination];
      vector<OligonucleotideFragment> vTmpFragmentList = OligoProcessUtils::GenerateOligonucleotideFragment(oligonucleotide, ionType, iCharge, varModCombination);
      vfragmentList.insert(vfragmentList.end(), vTmpFragmentList.begin(), vTmpFragmentList.end());
   }
   result.totalIons = vfragmentList.size();

   //generate theoretical spectrum and expeimental spectrum
   vector<double> theoreticalSpectrum;
   vector<double> experimentalSpectrum;

   bSuccessed = GenerateTheoreticalSpectrum(vfragmentList, expSpectrum.iArraySize, iCharge, theoreticalSpectrum);
   if (!bSuccessed)
       return false;
   bSuccessed = GenerateExperimentalSpectrum(expSpectrum, experimentalSpectrum);
   if (!bSuccessed)
       return false;

   // calculate Xcorr
   double dXcorr;
   bSuccessed = CalculateUtils::FastXcorrCalc(theoreticalSpectrum, experimentalSpectrum, Param::g_staticParams.iXcorrProcessingOffset, dXcorr);
   if (!bSuccessed)
       return false;
   result.xCorr = dXcorr;
   if (result.xCorr < CommonValues::dFloatZero)
      result.xCorr = 0;

   //cout match fragment
   vector<pair<int, int>> vIonSpectrumMatch = SearchUtils::FragmentMatchCount(expSpectrum, vfragmentList); //{peak, fragment}
   result.matchedIons = vIonSpectrumMatch.size();

   //calculate Sp bonus
   int iConsec = 0;
   for (int i = 0; i < vIonSpectrumMatch.size(); i++)
   {
      for (int j = i + 1; j < vIonSpectrumMatch.size(); j++)
      {
         if (vfragmentList[vIonSpectrumMatch[i].second].ionType == vfragmentList[vIonSpectrumMatch[j].second].ionType
         && vfragmentList[vIonSpectrumMatch[i].second].sSequence.size() == vfragmentList[vIonSpectrumMatch[j].second].sSequence.size() - 1)
            iConsec++;
      }
   }

   //calculate Sp
   double dMatchIntensity = 0.0;
   Spectrum spec = Spectrum(expSpectrum.spectrum);
   double dHighestIntensity = spec.at(expSpectrum.iHighestIonIndex).intensity;
   double dTotalIntensity = 0.0;
   for (auto& pair : vIonSpectrumMatch)
   {
      double dIntensity = spec.at(pair.first).intensity;
      dMatchIntensity += dIntensity;
   }

   result.dSp = dMatchIntensity/ expSpectrum.dTotalIntensity * ((double)(1+iConsec)/(double)(result.totalIons));
   if (result.dSp < CommonValues::dFloatZero)
      result.dSp = 0;
   return true;
}

bool SearchUtils::GenerateTheoreticalSpectrum(vector<OligonucleotideFragment> fragmentList, int iArraySize, int iMaxFragmentCharge, vector<double>& spectrum)
{
   spectrum.clear();
   spectrum.assign(iArraySize, 0.0);

   for (OligonucleotideFragment& fragment: fragmentList)
   {
      if (fragment.iCharge != 0)
      {
         double dMZ = MassUtils::MassToMZ(fragment.dMass, fragment.iCharge);
         int iBinIndex = BIN(dMZ);
         if (iBinIndex <= spectrum.size() - 1)
            spectrum[iBinIndex] += 1.0;
      }
      else
      {
         for (int charge = 1; charge <= abs(iMaxFragmentCharge); charge++)
         {
            double dMass = fragment.dMass;
            double dMZ = 0.0;
            if (iMaxFragmentCharge > 0)
               dMZ = MassUtils::MassToMZ(dMass, charge);
            else
               dMZ = MassUtils::MassToMZ(dMass, -1 * charge);
            int iBinIndex = BIN(dMZ);
            if (iBinIndex <= spectrum.size() - 1)
               spectrum[iBinIndex] += 1.0;
         }
      }
   }
   return true;
}

bool SearchUtils::GenerateExperimentalSpectrum(const ExpSpectrum& expSpectrum, vector<double>& vSpectrumSignal)
{
    vSpectrumSignal.clear();
    vSpectrumSignal.assign(expSpectrum.iArraySize, 0.0);
   Spectrum spectrum = Spectrum(expSpectrum.spectrum);
   double dMaxIntensity = 0.0;
   for (int i = 0; i < spectrum.size(); i++)
   {
      Peak_T peak = spectrum.at(i);

      if(peak.intensity > dMaxIntensity)
         dMaxIntensity = peak.intensity;

      double dMZ = peak.mz;
      int iBinIndex = BIN(dMZ);
      if (iBinIndex < vSpectrumSignal.size())
          vSpectrumSignal[iBinIndex] += peak.intensity;
   }

   for (int i = 0; i < vSpectrumSignal.size(); i++)
   {
       vSpectrumSignal[i] /= dMaxIntensity;
   }
   
   return true;
}

vector<pair<int, int>> SearchUtils::FragmentMatchCount(const ExpSpectrum& expSpectrum, const vector<OligonucleotideFragment>& vFragmentList)
{
   int iMatchCount = 0;
   vector<pair<int, int>> vMatchIndex;
   Spectrum spectrum = Spectrum(expSpectrum.spectrum);
   for (int i = 0; i < spectrum.size(); i++)
   {
      Peak_T peak = spectrum.at(i);
      if (peak.intensity < CommonValues::dFloatZero)
         continue;
      double dExpMZ = peak.mz;

      for (int ii = 0; ii < vFragmentList.size(); ii++)
      {
         const OligonucleotideFragment& fragment = vFragmentList[ii];
         double dMass = fragment.dMass;
         int iCharge = fragment.iCharge;

         if (Param::g_staticParams.tolerances.iFragmentToleranceUnits == 0) // amu
         {
            if (Param::g_staticParams.tolerances.iFragmentToleranceType == 1)  // precursor m/z tolerance
            {
               double dTheMZ = MassUtils::MassToMZ(dMass, iCharge);
               if (DoubleRange::inRange(Param::g_staticParams.tolerances.fragmentTolerance, dExpMZ - dTheMZ))
               {
                  iMatchCount ++;
                  vMatchIndex.push_back({i, ii});
                  break;
               }
            }
            else
            {
               double dExpMass = MassUtils::MZToMass(dExpMZ, iCharge);
               if (DoubleRange::inRange(Param::g_staticParams.tolerances.fragmentTolerance, dExpMass - dMass))
               {
                  iMatchCount ++;
                  vMatchIndex.push_back({i, ii});
                  break;
               }
            }
         }
         else if (Param::g_staticParams.tolerances.iFragmentToleranceUnits == 1) // mmu
         {
            DoubleRange tolerance = Param::g_staticParams.tolerances.fragmentTolerance;
            tolerance *= 0.001;
            if (Param::g_staticParams.tolerances.iFragmentToleranceType == 1)  // precursor m/z tolerance
            {
               double dTheMZ = MassUtils::MassToMZ(dMass, iCharge);
               if (DoubleRange::inRange(tolerance, dExpMZ - dTheMZ))
               {
                  iMatchCount ++;
                  vMatchIndex.push_back({i, ii});
                  break;
               }
            }
            else
            {
               double dExpMass = MassUtils::MZToMass(dExpMZ, iCharge);
               if (DoubleRange::inRange(tolerance, dExpMass - dMass))
               {
                  iMatchCount ++;
                  vMatchIndex.push_back({i, ii});
                  break;
               }
            }
         }
         else // ppm
         {
            double dTheMZ = MassUtils::MassToMZ(dMass, iCharge);

            DoubleRange tolerance = Param::g_staticParams.tolerances.fragmentTolerance;
            tolerance *= ((double)dTheMZ / (double)1e6);
            if (DoubleRange::inRange(tolerance, dExpMZ - dTheMZ))
            {
               iMatchCount ++;
               vMatchIndex.push_back({i, ii});
               break;
            }
         }
      }
   }
   return vMatchIndex;
}


