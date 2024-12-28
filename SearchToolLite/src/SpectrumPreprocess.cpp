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


#include "SpectrumPreprocess.h"
#include "MassUtils.h"


MSReader SpectrumPreprocessUtils::msreader;

SpectrumPreprocessUtils::SpectrumPreprocessUtils()  { }

SpectrumPreprocessUtils::~SpectrumPreprocessUtils() { }

MSReader& SpectrumPreprocessUtils::getMSReader()
{
    return SpectrumPreprocessUtils::msreader;
}

bool SpectrumPreprocessUtils::LoadAllSpectrum(vector<InputFile>& inputFiles, vector<std::pair<InputFile, vector<CometSpectrum>>>& vInputFileCometSpectrumPair)
{
    if (vInputFileCometSpectrumPair.size() > 0)
        vInputFileCometSpectrumPair.clear();

    for (InputFile inputFile: inputFiles)
    {
        vector<CometSpectrum> vCometSpectrumList;

        if(! SpectrumPreprocessUtils::LoadSpectrum(inputFile, vCometSpectrumList))
            return false;
        

        if (vCometSpectrumList.size() <= 0)
            return false;

        InputFile file = InputFile(inputFile);
        std::pair<InputFile, vector<CometSpectrum>> pair(file, vCometSpectrumList);
        vInputFileCometSpectrumPair.push_back(pair);
    }

    return true;
}

bool SpectrumPreprocessUtils::LoadSpectrum(InputFile& inputFile, vector<CometSpectrum>& vCometSpectrumList)
{
    int iFileLastScan = -1;                         // The actual last scan in the file.
    int iBreakMaxNum = 0;
    Spectrum mstSpectrum;                              // For holding spectrum.
    CometSpectrum cometSpectrum;

    SpectrumPreprocessUtils::msreader.readFile(inputFile.sFileName.c_str(), mstSpectrum, 0);       //first time read

    iFileLastScan = SpectrumPreprocessUtils::msreader.getLastScan();

    if ((iFileLastScan != -1) && (iFileLastScan < inputFile.iFirstScan))                        //invalid scan range
        return false;
    
    if (mstSpectrum.getScanNumber() < inputFile.iFirstScan)                                     //skip to the first scan
    {
        SpectrumPreprocessUtils::msreader.readFile(inputFile.sFileName.c_str(), mstSpectrum, inputFile.iFirstScan); 

        // iScanNumber will equal 0 if iFirstScan is not the right scan level so need to keep reading the next scan until we get a non-zero scan number
        while (mstSpectrum.getScanNumber() == 0 && inputFile.iFirstScan < inputFile.iLastScan)
        {
            inputFile.iFirstScan++;
            SpectrumPreprocessUtils::msreader.readFile(inputFile.sFileName.c_str(), mstSpectrum, inputFile.iFirstScan); 
        }

        if (mstSpectrum.getScanNumber() <= 0 || inputFile.iFirstScan > inputFile.iLastScan)
            return false;
    }

    while(true)
    {
        if (mstSpectrum.getScanNumber() != 0)   
        {
            // iFirstScan and iLastScan are both 0 unless scan range is specified. If -L last scan specified, iFirstScan is set to 1
            // However, if -F is specified, iLastscan is still 0. If scan range is specified, need to enforce here.
            if (inputFile.iFirstScan != 0 && inputFile.iLastScan != 0 && !(inputFile.iFirstScan <= mstSpectrum.getScanNumber() && mstSpectrum.getScanNumber() <= inputFile.iLastScan))
                continue;
            if (inputFile.iFirstScan != 0 && inputFile.iLastScan == 0 && mstSpectrum.getScanNumber() < inputFile.iFirstScan)
                continue;

            // Clear out m/z range if clear_mz_range parameter is specified. Accomplish this by setting corresponding intensity to 0
            if (Param::g_staticParams.options.clearMzRange.dEnd > 0.0 && Param::g_staticParams.options.clearMzRange.dStart <= Param::g_staticParams.options.clearMzRange.dEnd)
            {
                mstSpectrum.sortMZ();
                for (int i = 0; i < mstSpectrum.size(); i++)
                {
                    if (mstSpectrum.at(i).mz > Param::g_staticParams.options.clearMzRange.dEnd)
                        break;

                    if (Param::g_staticParams.options.clearMzRange.dStart <= mstSpectrum.at(i).mz && mstSpectrum.at(i).mz <= Param::g_staticParams.options.clearMzRange.dEnd)
                        mstSpectrum.at(i).intensity = 0.0;
                        
                }
            }

            if (inputFile.iAnalysisType == InputFile::AnalysisType_SpecificScanRange && inputFile.iLastScan > 0 && mstSpectrum.getScanNumber() > inputFile.iLastScan) // if specific scan then consider to break
                break;

            if (!SpectrumPreprocessUtils::PreprocessSpectrum(mstSpectrum, vCometSpectrumList))
                break;
            // cometSpectrum.spectrum = mstSpectrum;
            // spectrumList.push_back(cometSpectrum);
        }
        else
        {
            // What happens here when iScanNumber == 0 and it is an mzXML file?
            // Best way to deal with this is to keep trying to read but track each
            // attempt and break when the count goes past the mzXML's last scan.
            if (++iBreakMaxNum > iFileLastScan)
                break;
        }      

        SpectrumPreprocessUtils::msreader.readFile(nullptr, mstSpectrum);      
    }

    if (vCometSpectrumList.size() <= 0)
        return false;


    return true;
}


//根据Spectrum的不同前体离子电荷构建CometSpectrum
bool SpectrumPreprocessUtils::PreprocessSpectrum(Spectrum& spec, vector<CometSpectrum>& cometSpectrumList)
{
    // To run a search, all that's needed is MH+ and Z. So need to generate
    // all combinations of these for each spectrum, whether there's a known
    // Z for each precursor or if Comet has to guess the 1+ or 2+/3+ charges.
    for (int i = 0 ; i < spec.sizeMZ(); ++i)  // walk through all precursor m/z's; usually just one
    {
      double dMZ = 0.0;                      // m/z to use for analysis
      vector<int> vChargeStates;
      int iSpectrumCharge = 0;

      if ((dMZ = spec.getMonoMZ(i)) == 0)
         dMZ = spec.getMZ(i);

      if (spec.sizeMZ() <= spec.sizeZ())
         iSpectrumCharge = spec.atZ(i).z;
      else
         iSpectrumCharge = 0;

      // 1.  Have spectrum charge from file.  It may be 0.
      // 2.  If the precursor_charge range is set and override_charge is set, then do something look into charge range.
      // 3.  Else just use the spectrum charge (and possibly 1 or 2/3 rule)
      // for simple this vision use spectrum charge

      if (iSpectrumCharge != 0) // use charge from spectrum file
      {
         vector<int> vTmpChargeStates;
         vTmpChargeStates.push_back(iSpectrumCharge);

         // add in any other charge states for the single precursor m/z
         if (spec.sizeMZ() == 1 && spec.sizeMZ() < spec.sizeZ())
         {
               for (int ii = 1 ; ii < spec.sizeZ(); ++ii)
                  vTmpChargeStates.push_back(spec.atZ(ii).z);
         }

         for (auto& charge : vTmpChargeStates)
         {
            bool bPositive = Param::g_staticParams.options.bPositive;
            if ((bPositive && charge > 0) || (!bPositive && charge < 0))
               vChargeStates.push_back(charge);
         }
      }
      else
      {
         double dSumBelow = 0.0;
         double dSumTotal = 0.0;

         for (int i=0; i<spec.size(); ++i)
         {
            dSumTotal += spec.at(i).intensity;

            if (spec.at(i).mz < spec.getMZ())
               dSumBelow += spec.at(i).intensity;
         }

         if (isEqual(dSumTotal, 0.0) || ((dSumBelow/dSumTotal) > 0.95))
         {
            if (Param::g_staticParams.options.bPositive)
               vChargeStates.push_back(1);
            else
               vChargeStates.push_back(-1);
         }
         else
         {
            if (Param::g_staticParams.options.bPositive)
            {
               vChargeStates.push_back(2);
               vChargeStates.push_back(3);
            }
            else
            {
               vChargeStates.push_back(-2);
               vChargeStates.push_back(-3);
            }
         }
      }

      // now analyze all possible precursor charges for this spectrum
      for (vector<int>::iterator iter = vChargeStates.begin(); iter != vChargeStates.end(); ++iter)
      {
         int iPrecursorCharge = *iter;
         double dMass = MassUtils::MZToMass(dMZ, iPrecursorCharge);

         if (IntRange::inRange(Param::g_staticParams.options.precursorChargeRange, abs(iPrecursorCharge))
            && DoubleRange::inRange(Param::g_staticParams.options.oligonucleotideMassRange, dMass))
         {
            CometSpectrum cometSpectrum;
            cometSpectrum.spectrum = Spectrum(spec);
            cometSpectrum.dExperimentalPrecursorMass = dMass;
            cometSpectrum.iPrecursorChargeState = iPrecursorCharge;

            double dTotalIntensity = 0.0;
            double dHighestIntensity = 0;
            for (int j = 0; j < cometSpectrum.spectrum.size(); j++)
            {
               Peak_T peak = cometSpectrum.spectrum.at(j);
               if (peak.intensity > dHighestIntensity)
               {
                  dHighestIntensity = peak.intensity;
                  cometSpectrum.iHighestIonIndex = j;
               }
               dTotalIntensity += peak.intensity;
            }
            cometSpectrum.dTotalIntensity = dTotalIntensity;

            SpectrumPreprocessUtils::CalculateMassTolerance(cometSpectrum);

            cometSpectrum.iArraySize = (int)(cometSpectrum.expPrecursorMassSearchTolerance.dEnd * Param::g_staticParams.dInverseBinWidth);
            
            cometSpectrumList.push_back(cometSpectrum);
         }
      }
   }

   return true;
}

bool SpectrumPreprocessUtils::CalculateMassTolerance(CometSpectrum& cometSpectrum)
{
   if (Param::g_staticParams.tolerances.iPrecursorToleranceUnits == 0) // amu
   {
      cometSpectrum.expPrecursorMassTolerance = Param::g_staticParams.tolerances.precursorTolerance;
      if (Param::g_staticParams.tolerances.iPrecursorToleranceType == 1)  // precursor m/z tolerance
      {
         cometSpectrum.expPrecursorMassTolerance *= abs(cometSpectrum.iPrecursorChargeState);
      }
   }
   else if (Param::g_staticParams.tolerances.iPrecursorToleranceUnits == 1) // mmu
   {
      cometSpectrum.expPrecursorMassTolerance = Param::g_staticParams.tolerances.precursorTolerance;
      cometSpectrum.expPrecursorMassTolerance *= 0.001;
      if (Param::g_staticParams.tolerances.iPrecursorToleranceType == 1)  // precursor m/z tolerance
      {
         cometSpectrum.expPrecursorMassTolerance *= abs(cometSpectrum.iPrecursorChargeState);
      }
   }
   else // ppm
   {
      double dMZ = MassUtils::MassToMZ(cometSpectrum.dExperimentalPrecursorMass, cometSpectrum.iPrecursorChargeState);

      // calculate lower/upper ppm bounds in m/z
      double dMZLower = dMZ + (dMZ * Param::g_staticParams.tolerances.precursorTolerance.dStart / 1E6); // dInputToleranceMinus typically has negative sign
      double dMZUpper = dMZ + (dMZ * Param::g_staticParams.tolerances.precursorTolerance.dEnd  / 1E6);

      // convert m/z bounds to neutral mass then add a proton as Comet uses protonated ranges
      double dMassLower = MassUtils::MZToMass(dMZLower, cometSpectrum.iPrecursorChargeState);
      double dMassUpper = MassUtils::MZToMass(dMZUpper, cometSpectrum.iPrecursorChargeState);

      // these are now the mass difference from exp mass (experimental M mass) to lower/upper bounds
      cometSpectrum.expPrecursorMassTolerance.dStart  = dMassLower - cometSpectrum.dExperimentalPrecursorMass; // usually negative
      cometSpectrum.expPrecursorMassTolerance.dEnd = dMassUpper - cometSpectrum.dExperimentalPrecursorMass; // usually positive
   }

   double dC13Diff = (ElementMassUtils::Carbon_13_Mono - ElementMassUtils::Carbon_Mono);
   vector<int> vC13DiffNums;
   vC13DiffNums.push_back(0);
   if (Param::g_staticParams.tolerances.iIsotopeSearchType == 0) //search 0 isotope windows
   {
      cometSpectrum.expPrecursorMassSearchTolerance.dStart = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dStart;
      cometSpectrum.expPrecursorMassSearchTolerance.dEnd = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dEnd;
   }
   else if (Param::g_staticParams.tolerances.iIsotopeSearchType == 1) // search 0, +1 isotope windows
   {
      cometSpectrum.expPrecursorMassSearchTolerance.dStart = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dStart - dC13Diff;
      cometSpectrum.expPrecursorMassSearchTolerance.dEnd = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dEnd;
      vC13DiffNums.push_back(-1);
   }
   else if (Param::g_staticParams.tolerances.iIsotopeSearchType == 2) // search 0, +1, +2 isotope windows
   {
      cometSpectrum.expPrecursorMassSearchTolerance.dStart = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dStart - 2.0 * dC13Diff;
      cometSpectrum.expPrecursorMassSearchTolerance.dEnd = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dEnd;
      vC13DiffNums.push_back(-1);
      vC13DiffNums.push_back(-2);
   }
   else if (Param::g_staticParams.tolerances.iIsotopeSearchType == 3) // search 0, +1, +2, +3 isotope windows
   {
      cometSpectrum.expPrecursorMassSearchTolerance.dStart = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dStart - 3.0 * dC13Diff;
      cometSpectrum.expPrecursorMassSearchTolerance.dEnd = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dEnd;
      vC13DiffNums.push_back(-1);
      vC13DiffNums.push_back(-2);
      vC13DiffNums.push_back(-3);
   }
   else if (Param::g_staticParams.tolerances.iIsotopeSearchType == 4) // search -1, 0, +1, +2, +3 isotope windows
   {
      cometSpectrum.expPrecursorMassSearchTolerance.dStart = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dStart - 3.0 * dC13Diff;
      cometSpectrum.expPrecursorMassSearchTolerance.dEnd = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dEnd + 1.0 * dC13Diff;
      vC13DiffNums.push_back(-1);
      vC13DiffNums.push_back(-2);
      vC13DiffNums.push_back(-3);
      vC13DiffNums.push_back(1);
   }
   else if (Param::g_staticParams.tolerances.iIsotopeSearchType == 5) // search -1, 0, +1 isotope windows
   {
      cometSpectrum.expPrecursorMassSearchTolerance.dStart = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dStart - dC13Diff;
      cometSpectrum.expPrecursorMassSearchTolerance.dEnd = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dEnd + dC13Diff;
      vC13DiffNums.push_back(-1);
      vC13DiffNums.push_back(1);
   }
   else if (Param::g_staticParams.tolerances.iIsotopeSearchType == 6) // search -3, -2, -1, 0, +1, +2, +3 isotope windows
   {
      cometSpectrum.expPrecursorMassSearchTolerance.dStart = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dStart - 3.0 * dC13Diff;
      cometSpectrum.expPrecursorMassSearchTolerance.dEnd = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dEnd + 3.0 * dC13Diff;
      vC13DiffNums.push_back(1);
      vC13DiffNums.push_back(2);
      vC13DiffNums.push_back(3);
      vC13DiffNums.push_back(-1);
      vC13DiffNums.push_back(-2);
      vC13DiffNums.push_back(-3);
   }
   else if (Param::g_staticParams.tolerances.iIsotopeSearchType == 7) // search -8, -4, 0, 4, 8 windows
   {
      cometSpectrum.expPrecursorMassSearchTolerance.dStart = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dStart - 8.0 * dC13Diff;
      cometSpectrum.expPrecursorMassSearchTolerance.dEnd = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dEnd + 8.0 * dC13Diff;
      vC13DiffNums.push_back(-8);
      vC13DiffNums.push_back(-4);
      vC13DiffNums.push_back(4);
      vC13DiffNums.push_back(8);
   }
   else  // Should not get here.
   {
      char szErrorMsg[256];
      sprintf(szErrorMsg,  " Error - iIsotopeError=%d\n",  Param::g_staticParams.tolerances.iIsotopeSearchType);
      logerr(szErrorMsg);
      cometSpectrum.expPrecursorMassSearchTolerance.dStart = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dStart;
      cometSpectrum.expPrecursorMassSearchTolerance.dEnd = cometSpectrum.dExperimentalPrecursorMass + cometSpectrum.expPrecursorMassTolerance.dEnd;
   }

   // now set Low/High to mass range around ExpMass
   cometSpectrum.expPrecursorMassTolerance += cometSpectrum.dExperimentalPrecursorMass;

   for (int& num : vC13DiffNums)
   {
      DoubleRange range = cometSpectrum.expPrecursorMassTolerance;
      range += num * dC13Diff;
      cometSpectrum.vExpPrecursorMassSearchRange.push_back(range);
   }

   return true;
}



