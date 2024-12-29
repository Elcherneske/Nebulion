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

SpectrumPreprocessUtils::SpectrumPreprocessUtils()  
{

}

SpectrumPreprocessUtils::~SpectrumPreprocessUtils() 
{ 
       
}

MSReader& SpectrumPreprocessUtils::getMSReader()
{
    return SpectrumPreprocessUtils::msreader;
}

bool SpectrumPreprocessUtils::LoadSpectrum(InputFile& inputFile, vector<ExpSpectrum>& vSpectrumList)
{
    int iFileLastScan = -1;                         // The actual last scan in the file.
    int iBreakMaxNum = 0;
    Spectrum mstSpectrum;                              // For holding spectrum.
    ExpSpectrum expSpectrum;

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

            if (!SpectrumPreprocessUtils::PreprocessSpectrum(mstSpectrum, vSpectrumList))
                break;
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

    if (vSpectrumList.size() <= 0)
        return false;


    return true;
}


//根据Spectrum的不同前体离子电荷构建expSpectrum
bool SpectrumPreprocessUtils::PreprocessSpectrum(Spectrum& spec, vector<ExpSpectrum>& vSpectrumList)
{
    // To run a search, all that's needed is MZ and Z. So need to generate all combinations of these for each spectrum, whether there's a known
    // Z for each precursor or add all posible Z to let the MZ in range
    for (int mz_index = 0 ; mz_index < spec.sizeMZ(); ++mz_index)  // walk through all precursor m/z's; usually just one
    {
        double dMZ = 0.0;                                   // m/z to use for analysis
        vector<int> vChargeStates;
        int iSpectrumCharge = 0;

        if ((dMZ = spec.getMonoMZ(mz_index)) == 0)
            dMZ = spec.getMZ(mz_index);

        if ( mz_index < spec.sizeZ() )
            iSpectrumCharge = spec.atZ(mz_index).z;
        else
            iSpectrumCharge = 0;

      if (iSpectrumCharge != 0) // use charge from spectrum file
      {
         vector<int> vTmpChargeStates;
         vTmpChargeStates.push_back(iSpectrumCharge);

         // add in any other charge states for the single precursor m/z
         if (spec.sizeMZ() == 1 && spec.sizeMZ() < spec.sizeZ())
         {
            for (int i = 1; i < spec.sizeZ(); ++i)
                vTmpChargeStates.push_back(spec.atZ(i).z);
         }

         for (auto& charge : vTmpChargeStates)
         {
            bool bPositive = Param::g_staticParams.options.bPositive;
            if ((bPositive && charge > 0) || (!bPositive && charge < 0))
                vChargeStates.push_back(charge);
            else
                vChargeStates.push_back(-1 * charge);
         }
      }
      else
      {
         double dSumBelow = 0.0;
         double dSumTotal = 0.0;
         for (int i=0; i<spec.size(); ++i)
         {
            dSumTotal += spec.at(i).intensity;
            if (spec.at(i).mz < dMZ)
               dSumBelow += spec.at(i).intensity;
         }

         if ((dSumTotal < CommonValues::dFloatZero) || ((dSumBelow/dSumTotal) > 0.95))
         {
            if (Param::g_staticParams.options.bPositive)
               vChargeStates.push_back(1);
            else
               vChargeStates.push_back(-1);
         }
         else
         {
             IntRange chargeRange = Param::g_staticParams.options.precursorChargeRange;
            if (Param::g_staticParams.options.bPositive)
            {
                for (int charge = chargeRange.iStart; charge <= chargeRange.iEnd; charge++)
                {
                    vChargeStates.push_back(charge);
                }
            }
            else
            {
                for (int charge = chargeRange.iStart; charge <= chargeRange.iEnd; charge++)
                {
                    vChargeStates.push_back(-1 * charge);
                }
            }
         }
      }

      // now analyze all possible precursor charges for this spectrum
      for (vector<int>::iterator iter = vChargeStates.begin(); iter != vChargeStates.end(); ++iter)
      {
         int iPrecursorCharge = *iter;
         double dMass = MassUtils::MZToMass(dMZ, iPrecursorCharge);

         if (Param::g_staticParams.options.precursorChargeRange.inRange(abs(iPrecursorCharge))
            && Param::g_staticParams.options.oligonucleotideMassRange.inRange(dMass))
         {
             ExpSpectrum expSpectrum;
             expSpectrum.spectrum = Spectrum(spec);
             expSpectrum.dExpPreMass = dMass;
             expSpectrum.iPreCharge = iPrecursorCharge;

            double dTotalIntensity = 0.0;
            double dHighestIntensity = 0;
            for (int j = 0; j < expSpectrum.spectrum.size(); j++)
            {
               Peak_T peak = expSpectrum.spectrum.at(j);
               if (peak.intensity > dHighestIntensity)
               {
                  dHighestIntensity = peak.intensity;
                  expSpectrum.iHighestIonIndex = j;
               }
               dTotalIntensity += peak.intensity;
            }
            expSpectrum.dTotalIntensity = dTotalIntensity;

            SpectrumPreprocessUtils::CalculateMassTolerance(expSpectrum);

            expSpectrum.iArraySize = (int)(expSpectrum.preMassSearchRange.dEnd * Param::g_staticParams.dInverseBinWidth);
            
            vSpectrumList.push_back(expSpectrum);
         }
      }
   }

   return true;
}

bool SpectrumPreprocessUtils::CalculateMassTolerance(ExpSpectrum& expSpectrum)
{
   if (Param::g_staticParams.tolerances.iPrecursorToleranceUnits == 0) // amu
   {
       expSpectrum.preBasicTolerance = Param::g_staticParams.tolerances.precursorTolerance;
      if (Param::g_staticParams.tolerances.iPrecursorToleranceType == 1)  // precursor m/z tolerance
      {
          expSpectrum.preBasicTolerance *= abs(expSpectrum.iPreCharge);
      }
   }
   else if (Param::g_staticParams.tolerances.iPrecursorToleranceUnits == 1) // mmu
   {
       expSpectrum.preBasicTolerance = Param::g_staticParams.tolerances.precursorTolerance;
       expSpectrum.preBasicTolerance *= 0.001;
      if (Param::g_staticParams.tolerances.iPrecursorToleranceType == 1)  // precursor m/z tolerance
      {
          expSpectrum.preBasicTolerance *= abs(expSpectrum.iPreCharge);
      }
   }
   else // ppm
   {
      double dMZ = MassUtils::MassToMZ(expSpectrum.dExpPreMass, expSpectrum.iPreCharge);

      // calculate lower/upper ppm bounds in m/z
      double dMZLower = dMZ + (dMZ * Param::g_staticParams.tolerances.precursorTolerance.dStart / 1E6); // dInputToleranceMinus typically has negative sign
      double dMZUpper = dMZ + (dMZ * Param::g_staticParams.tolerances.precursorTolerance.dEnd  / 1E6);

      double dMassLower = MassUtils::MZToMass(dMZLower, expSpectrum.iPreCharge);
      double dMassUpper = MassUtils::MZToMass(dMZUpper, expSpectrum.iPreCharge);

      // these are now the mass difference from exp mass (experimental M mass) to lower/upper bounds
      expSpectrum.preBasicTolerance.dStart  = dMassLower - expSpectrum.dExpPreMass; // usually negative
      expSpectrum.preBasicTolerance.dEnd = dMassUpper - expSpectrum.dExpPreMass; // usually positive
   }

   double dC13Diff = (ElementMassUtils::Carbon_13_Mono - ElementMassUtils::Carbon_Mono);
   double dN15Diff = (ElementMassUtils::Nitrogen_15_Mono - ElementMassUtils::Nitrogen_Mono);
   vector<int> vDiffNums;
   vDiffNums.push_back(0);
   if (Param::g_staticParams.tolerances.iIsotopeSearchType == 0) //search 0 isotope windows
   {

   }
   else if (Param::g_staticParams.tolerances.iIsotopeSearchType == 1) // search 0, -1 isotope windows
   {
      vDiffNums.push_back(-1);
   }
   else if (Param::g_staticParams.tolerances.iIsotopeSearchType == 2) // search 0, -1, -2 isotope windows
   {
       vDiffNums.push_back(-1);
       vDiffNums.push_back(-2);
   }
   else if (Param::g_staticParams.tolerances.iIsotopeSearchType == 3) // search 0, -1, -2, -3 isotope windows
   {
       vDiffNums.push_back(-1);
       vDiffNums.push_back(-2);
       vDiffNums.push_back(-3);
   }
   else if (Param::g_staticParams.tolerances.iIsotopeSearchType == 4) // search +1, 0, -1, -2, -3 isotope windows
   {
      vDiffNums.push_back(-1);
      vDiffNums.push_back(-2);
      vDiffNums.push_back(-3);
      vDiffNums.push_back(1);
   }
   else if (Param::g_staticParams.tolerances.iIsotopeSearchType == 5) // search +1, 0, -1 isotope windows
   {
       vDiffNums.push_back(-1);
       vDiffNums.push_back(1);
   }
   else if (Param::g_staticParams.tolerances.iIsotopeSearchType == 6) // search -3, -2, -1, 0, +1, +2, +3 isotope windows
   {
       vDiffNums.push_back(1);
       vDiffNums.push_back(2);
       vDiffNums.push_back(3);
       vDiffNums.push_back(-1);
       vDiffNums.push_back(-2);
       vDiffNums.push_back(-3);
   }
   else if (Param::g_staticParams.tolerances.iIsotopeSearchType == 7) // search -8, -4, 0, 4, 8 windows
   {
       vDiffNums.push_back(-8);
       vDiffNums.push_back(-4);
       vDiffNums.push_back(4);
       vDiffNums.push_back(8);
   }
   else  // Should not get here.
   {
       cout << "Error - iIsotopeError = " << Param::g_staticParams.tolerances.iIsotopeSearchType << endl;
   }

   // now set Low/High to mass range around ExpMass
   
   expSpectrum.preBasicTolerance += expSpectrum.dExpPreMass;

   for (int& num : vDiffNums)
   {
       DoubleRange range;
       range = expSpectrum.preBasicTolerance;
      range += num * dC13Diff;
      expSpectrum.vPreMassSearchRange.push_back(range);

      range = expSpectrum.preBasicTolerance;
      range += num * dN15Diff;
      expSpectrum.vPreMassSearchRange.push_back(range);
   }

   double dSearchMassStart = expSpectrum.preBasicTolerance.dStart;
   double dSearchMassEnd = expSpectrum.preBasicTolerance.dEnd;
   for (auto& range : expSpectrum.vPreMassSearchRange)
   {
       if (range.dStart < dSearchMassStart)
           dSearchMassStart = range.dStart;

       if (range.dEnd > dSearchMassEnd)
           dSearchMassEnd = range.dEnd;
   }

   expSpectrum.preMassSearchRange.dStart = dSearchMassStart;
   expSpectrum.preMassSearchRange.dEnd = dSearchMassEnd;

   return true;
}



