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

SpectrumPreprocessUtils::SpectrumPreprocessUtils()  
{

}

SpectrumPreprocessUtils::~SpectrumPreprocessUtils() 
{ 
       
}

MSReader& SpectrumPreprocessUtils::GetMSReader()
{
    return this->msreader;
}

bool SpectrumPreprocessUtils::LoadSpectrum(InputFile& inputFile, vector<ExpSpectrum>& vSpectrumList)
{
    int iFileLastScan = -1;                         // The actual last scan in the file.
    int iBreakMaxNum = 0;
    Spectrum mstSpectrum;                              // For holding spectrum.
    ExpSpectrum expSpectrum;

    this->msreader.readFile(inputFile.sFileName.c_str(), mstSpectrum, 0);       //first time read
    iFileLastScan = this->msreader.getLastScan();
    if ((iFileLastScan != -1) && (iFileLastScan < inputFile.iFirstScan))                        //invalid scan range
        return false;
    
    if (mstSpectrum.getScanNumber() < inputFile.iFirstScan)                                     //skip to the first scan
    {
        this->msreader.readFile(inputFile.sFileName.c_str(), mstSpectrum, inputFile.iFirstScan); 
        // iScanNumber will equal 0 if iFirstScan is not the right scan level so need to keep reading the next scan until we get a non-zero scan number
        while (mstSpectrum.getScanNumber() == 0 && inputFile.iFirstScan < inputFile.iLastScan)
        {
            inputFile.iFirstScan++;
            this->msreader.readFile(inputFile.sFileName.c_str(), mstSpectrum, inputFile.iFirstScan); 
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

            if (inputFile.iAnalysisType == InputFile::AnalysisType_SpecificScanRange && inputFile.iLastScan > 0 && mstSpectrum.getScanNumber() > inputFile.iLastScan) // if specific scan then consider to break
                break;

            if (!this->PreprocessSpectrum(mstSpectrum, vSpectrumList))
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

        this->msreader.readFile(nullptr, mstSpectrum);      
    }

    if (vSpectrumList.size() <= 0)
        return false;

    return true;
}


//根据Spectrum的不同前体离子电荷构建expSpectrum
bool SpectrumPreprocessUtils::PreprocessSpectrum(Spectrum& spec, vector<ExpSpectrum>& vSpectrumList)
{
    if (spec.sizeMZ() == 0)
        return false;
    
    double dMZ = spec.getMonoMZ(0);
    if (dMZ == 0)
        dMZ = spec.getMZ(0);
    
    int iSpectrumCharge;
    if (spec.sizeZ() > 0)
        iSpectrumCharge = spec.atZ(0).z;
    else
        iSpectrumCharge = 0;

    // generate all possible charge states
    vector<int> vChargeStates;
    if (iSpectrumCharge != 0) // use charge from spectrum file
    {
        vector<int> vTmpChargeStates;
        vTmpChargeStates.push_back(iSpectrumCharge);
        // add in any other charge states for the single precursor m/z
        if (spec.sizeZ() > 1)
        {
            for (int i = 1; i < spec.sizeZ(); ++i)
                vTmpChargeStates.push_back(spec.atZ(i).z);
        }
        for (auto& charge : vTmpChargeStates)
        {
            bool bPositive = Param::getInstance().options.bPositiveCharge;
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
            Param::getInstance().options.bPositiveCharge ? vChargeStates.push_back(1) : vChargeStates.push_back(-1);
        else
        {
            IntRange chargeRange = Param::getInstance().options.precursorChargeRange;
            if (Param::getInstance().options.bPositiveCharge)
                for (int charge = chargeRange.iStart; charge <= chargeRange.iEnd; charge++)
                    vChargeStates.push_back(charge);
            else
                for (int charge = chargeRange.iStart; charge <= chargeRange.iEnd; charge++)
                    vChargeStates.push_back(-1 * charge);
        }
    }

    // now analyze all possible precursor charges for this spectrum
    for (int iPrecursorCharge : vChargeStates)
    {
        // double dMass = MassUtils::getInstance().MZToMass(dMZ, iPrecursorCharge, Param::getInstance().options.sAdductIon);
        if (Param::getInstance().options.precursorChargeRange.inRange(abs(iPrecursorCharge)) /*&& Param::getInstance().options.oligonucleotideMassRange.inRange(dMass)*/)
        {
            ExpSpectrum expSpectrum;
            expSpectrum.dExpMZ = dMZ;
            expSpectrum.iPreCharge = iPrecursorCharge;

            // calculate highest intensity
            double dHighestIntensity = 0;
            for (int j = 0; j < spec.size(); j++)
            {
                Peak_T peak = spec.at(j);
                if (peak.intensity > dHighestIntensity)
                    dHighestIntensity = peak.intensity;
            }

            //generate vector<Peak> spectrum
            for (int j = 0; j < spec.size(); j++)
            {
                Peak_T peak = spec.at(j);
                if (peak.intensity < Param::getInstance().options.dMinimumIntensity || peak.intensity < dHighestIntensity * Param::getInstance().options.dPercentageBasePeak)
                    continue;
                if (Param::getInstance().options.bRemovePrecursorPeak && fabs(peak.mz - dMZ) < Param::getInstance().options.dRemovePrecursorTolerance)
                    continue;
                if (Param::getInstance().options.clearMZRange.inRange(peak.mz))
                    continue;
                
                expSpectrum.spectrum.push_back(ExpSpectrum::Peak(peak.mz, peak.intensity));
            }

            // calculate total intensity && highest intensity index
            double dTotalIntensity = 0;
            dHighestIntensity = 0;
            for (size_t j = 0; j < expSpectrum.spectrum.size(); j++)
            {
                dTotalIntensity += expSpectrum.spectrum[j].dIntensity;
                if (expSpectrum.spectrum[j].dIntensity > dHighestIntensity)
                    expSpectrum.iHighestIonIndex = j;
            }
            expSpectrum.dTotalIntensity = dTotalIntensity;

            this->CalculateMassTolerance(expSpectrum);

            expSpectrum.iArraySize = (int)(expSpectrum.preMassSearchRange.dEnd / Param::getInstance().options.dBinWidth);
            
            vSpectrumList.push_back(expSpectrum);
        }
    }

   return true;
}

bool SpectrumPreprocessUtils::CalculateMassTolerance(ExpSpectrum& expSpectrum)
{

    double dMZ = expSpectrum.dExpMZ;

    // calculate lower/upper ppm bounds in m/z
    DoubleRange mzTolerance = Param::getInstance().options.precursorTolerance;
    mzTolerance *= (dMZ/1E6);
    mzTolerance += dMZ;

    DoubleRange massTolerance;
    massTolerance.dStart = MassUtils::getInstance().MZToMass(mzTolerance.dStart, expSpectrum.iPreCharge, Param::getInstance().options.sAdductIon);
    massTolerance.dEnd = MassUtils::getInstance().MZToMass(mzTolerance.dEnd, expSpectrum.iPreCharge, Param::getInstance().options.sAdductIon);

   double dC13Diff = (ElementMassUtils::Carbon_13_Mono - ElementMassUtils::Carbon_Mono);
   double dN15Diff = (ElementMassUtils::Nitrogen_15_Mono - ElementMassUtils::Nitrogen_Mono);
   vector<int> vDiffNums;
   vDiffNums.push_back(0);
   if (Param::getInstance().options.iIsotopeError == 0) //search 0 isotope windows
   {

   }
   else if (Param::getInstance().options.iIsotopeError == 1) // search 0, -1 isotope windows
   {
      vDiffNums.push_back(-1);
   }
   else if (Param::getInstance().options.iIsotopeError == 2) // search 0, -1, -2 isotope windows
   {
       vDiffNums.push_back(-1);
       vDiffNums.push_back(-2);
   }
   else if (Param::getInstance().options.iIsotopeError == 3) // search 0, -1, -2, -3 isotope windows
   {
       vDiffNums.push_back(-1);
       vDiffNums.push_back(-2);
       vDiffNums.push_back(-3);
   }
   else if (Param::getInstance().options.iIsotopeError == 4) // search +1, 0, -1, -2, -3 isotope windows
   {
      vDiffNums.push_back(-1);
      vDiffNums.push_back(-2);
      vDiffNums.push_back(-3);
      vDiffNums.push_back(1);
   }
   else if (Param::getInstance().options.iIsotopeError == 5) // search +1, 0, -1 isotope windows
   {
       vDiffNums.push_back(-1);
       vDiffNums.push_back(1);
   }
   else if (Param::getInstance().options.iIsotopeError == 6) // search -3, -2, -1, 0, +1, +2, +3 isotope windows
   {
       vDiffNums.push_back(1);
       vDiffNums.push_back(2);
       vDiffNums.push_back(3);
       vDiffNums.push_back(-1);
       vDiffNums.push_back(-2);
       vDiffNums.push_back(-3);
   }
   else if (Param::getInstance().options.iIsotopeError == 7) // search -8, -4, 0, 4, 8 windows
   {
       vDiffNums.push_back(-8);
       vDiffNums.push_back(-4);
       vDiffNums.push_back(4);
       vDiffNums.push_back(8);
   }
   else  // Should not get here.
   {
       std::cerr << "Error - iIsotopeError = " << Param::getInstance().options.iIsotopeError << endl;
       throw std::runtime_error("Error - iIsotopeError = " + std::to_string(Param::getInstance().options.iIsotopeError));
   }

   // now set Low/High to mass range around ExpMass
   for (int& num : vDiffNums)
   {
        if (num == 0)
        {
            expSpectrum.vPreMassMatchRange.push_back(massTolerance);
            continue;
        }
        DoubleRange range;
        //C13
        range = massTolerance;
        range += num * dC13Diff;
        expSpectrum.vPreMassMatchRange.push_back(range);

        //N15
        range = massTolerance;
        range += num * dN15Diff;
        expSpectrum.vPreMassMatchRange.push_back(range);
   }

   if (expSpectrum.vPreMassMatchRange.size() <= 0) throw std::runtime_error("Error - No mass tolerance range found");

   //generate preMassSearchRange
   expSpectrum.preMassSearchRange = expSpectrum.vPreMassMatchRange[0];
   for (auto& range : expSpectrum.vPreMassMatchRange)
   {
       if (range.dStart < expSpectrum.preMassSearchRange.dStart)
           expSpectrum.preMassSearchRange.dStart = range.dStart;

       if (range.dEnd > expSpectrum.preMassSearchRange.dEnd)
           expSpectrum.preMassSearchRange.dEnd = range.dEnd;
   }

   return true;
}



