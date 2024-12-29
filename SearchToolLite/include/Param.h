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


#ifndef _PARAM_H_
#define _PARAM_H_


#include "Common.h"
#include "DataStructure.h"
#include "Threading.h"

struct Options 
{
   // charge
   IntRange precursorChargeRange;
   IntRange fragmentChargeRange;
   bool     bPositive;
   string   sAdductIon;

   //thread
   int iNumThreads;                    // 0=poll CPU else set # threads to spawn
   int iNumFragmentThreads;            // # threads used for fragment indexing

   //Oligonucleotide && Nucleic Acid
   IntRange oligonucleotideLengthRange;
   DoubleRange oligonucleotideMassRange;          // MH+ mass
   DoubleRange fragmentMassRange;                 // fragment index maximum fragment mass


   //spectrum
   IntRange    scanRange;
   DoubleRange clearMzRange;
   double dMinIntensity;                  // intensity cutoff for each peak
   double dMinPercentageIntensity;        // intensity cutoff for each peak as % of base peak
   int iSpectrumBatchSize;                // # of spectra to search at a time within the scan range
   int iWhichReadingFrame;
   int iMSLevel;                          // mzXML only
   int iMinPeaks;

   //static result
   double dMinimumXcorr;               // set the minimum xcorr to report (default is 1e-8)
   int bExplicitDeltaCn;               // if set to 1, do not use sequence similarity logic

   //Various Modification
   int iMaxVarModPerOligonucleoitde;   

   //others
   int iDecoySearch;                   // 0=no, 1=concatenated search, 2=separate decoy search
   int iDecoyKind;                     // 0=reverse, 1=shuffle
   

   Options() {
        precursorChargeRange = IntRange(0, CommonValues::iMaxOligonucleotideCharge);
        fragmentChargeRange  = IntRange(0, CommonValues::iMaxFragmentCharge);
        bPositive = false;
        sAdductIon = "H-"; //默认失去H离子，为负离子模式
        iNumThreads = 4;
        iNumFragmentThreads = 4;
        oligonucleotideLengthRange = IntRange(1, CommonValues::iMaxOligonucleotideLength);   
        oligonucleotideMassRange = DoubleRange(600.0, 8000.0);   
        fragmentMassRange = DoubleRange(100.0, 2000.0);  
        scanRange = IntRange(0, 0);              
        clearMzRange = DoubleRange(0.0, 0.0);     
        dMinIntensity = 0.0;
        dMinPercentageIntensity = 0.0;
        iSpectrumBatchSize = 0;
        iWhichReadingFrame = 0;
        iMSLevel = 2;
        iMinPeaks = 10;
        dMinimumXcorr = 1e-6;      
        bExplicitDeltaCn = 0;
        iMaxVarModPerOligonucleoitde = 3;
        iDecoySearch = 0;
        iDecoyKind = 0;
    }

   Options& operator=(const Options& other) {
      if (this != &other) {
         // Charge
         precursorChargeRange = other.precursorChargeRange;
         fragmentChargeRange = other.fragmentChargeRange;
         bPositive = other.bPositive;
         sAdductIon = other.sAdductIon;

         // Thread
         iNumThreads = other.iNumThreads;
         iNumFragmentThreads = other.iNumFragmentThreads;

         // Oligonucleotide && Nucleic Acid
         oligonucleotideLengthRange = other.oligonucleotideLengthRange;
         oligonucleotideMassRange = other.oligonucleotideMassRange;
         fragmentMassRange = other.fragmentMassRange;

         // Spectrum
         scanRange = other.scanRange;
         clearMzRange = other.clearMzRange;
         dMinIntensity = other.dMinIntensity;
         dMinPercentageIntensity = other.dMinPercentageIntensity;
         iSpectrumBatchSize = other.iSpectrumBatchSize;
         iWhichReadingFrame = other.iWhichReadingFrame;
         iMSLevel = other.iMSLevel;
         iMinPeaks = other.iMinPeaks;

         // Static result
         dMinimumXcorr = other.dMinimumXcorr;
         bExplicitDeltaCn = other.bExplicitDeltaCn;

         //various Modifications
         iMaxVarModPerOligonucleoitde = other.iMaxVarModPerOligonucleoitde;

         // Others
         iDecoySearch = other.iDecoySearch;
         iDecoyKind = other.iDecoyKind;
      }
      return *this;
   }
};


struct DBInfo //Fasta库信息
{
   string            sDatabase;                    //xxx.fasta.idx
   string            sFilename;                    //xxx.fasta
   int               iTotalNumNucleicAcid;
   unsigned long int uliTotNTCount;

   DBInfo()
   {
      sDatabase = "";
      sFilename = "";
      iTotalNumNucleicAcid = uliTotNTCount = 0;
   }

   DBInfo& operator=(const DBInfo& a)
   {
      sDatabase = a.sDatabase;
      sFilename = a.sFilename;
      iTotalNumNucleicAcid = a.iTotalNumNucleicAcid;
      uliTotNTCount = a.uliTotNTCount;
      return *this;
   }
};

// The minimum and maximum mass range of all oligonucleotide to consider and all fragment to consider, prepare to delete
struct MassChargeRange
{
   DoubleRange precursorMassRange;
   DoubleRange fragmentMassRange;
   IntRange precursorChargeRange;
   IntRange fragmentChargeRange;

   MassChargeRange()
   {
      precursorMassRange = DoubleRange(0.0, 0.0);
      fragmentMassRange = DoubleRange(0.0, 0.0);
      precursorChargeRange = IntRange(0, 0);
      fragmentChargeRange = IntRange(0, 0);
   }

   MassChargeRange& operator= (const MassChargeRange& other) {
      if (this != &other) {
         precursorMassRange = other.precursorMassRange;
         fragmentMassRange = other.fragmentMassRange;
         precursorChargeRange = other.precursorChargeRange;
         fragmentChargeRange = other.fragmentChargeRange;
      }
      return *this;
   }
   
};

struct ModificationInfo
{
   vector<Modification> vAdd5terminusOligonucleotide;
   vector<Modification> vAdd3terminusOligonucleotide;
   vector<Modification> vAdd5terminusNucleicAcid;
   vector<Modification> vAdd3terminusNucleicAcid;
   Modification pdStaticMods[CommonValues::iAsciiValueSize];
   string sStaticModificationFilename;
   string sVariableModificationFilename;
   string sModificationListFilename;

   ModificationInfo() 
   {
      for (int i = 0; i < CommonValues::iAsciiValueSize; ++i)
         pdStaticMods[i] = Modification();
   }

   ModificationInfo& operator=(const ModificationInfo& a)
   {
      vAdd5terminusOligonucleotide = a.vAdd5terminusOligonucleotide;
      vAdd3terminusOligonucleotide = a.vAdd3terminusOligonucleotide;
      vAdd5terminusNucleicAcid = a.vAdd5terminusNucleicAcid;
      vAdd3terminusNucleicAcid = a.vAdd3terminusNucleicAcid;
      sStaticModificationFilename = a.sStaticModificationFilename;
      sVariableModificationFilename = a.sVariableModificationFilename;
      sModificationListFilename = a.sModificationListFilename;

      for (int i = 0; i < CommonValues::iAsciiValueSize; ++i)
      {
         pdStaticMods[i] = a.pdStaticMods[i];
      }

      return *this;
   }
};

struct MassUtil
{
   int    bMonoMass;
   int    bDeoxidation;
   double dBackbone;
   double pdAtomMass[CommonValues::iAsciiValueSize];
   double pdNTMass[CommonValues::iAsciiValueSize];
   double pdBaseMass[CommonValues::iAsciiValueSize];
   double pdNTMassUser[CommonValues::iAsciiValueSize];       // user defined default amino acid masses

   MassUtil()
   {
      bMonoMass = 1;
      bDeoxidation = 0;
      dBackbone = 0.0;
      for (int i = 0; i < CommonValues::iAsciiValueSize; ++i)
      {
         pdAtomMass[i] = 999999.;
         pdBaseMass[i] = 999999.;
         pdNTMass[i] = 999999.;
         pdNTMassUser[i] = 0.0;
      }
   }

   MassUtil& operator=(const MassUtil& a)
   {
      bMonoMass = a.bMonoMass;
      bDeoxidation = a.bDeoxidation;
      dBackbone = a.dBackbone;
      for (int i = 0; i < CommonValues::iAsciiValueSize; ++i)
      {
         pdAtomMass[i] = a.pdAtomMass[i];
         pdNTMass[i] = a.pdNTMass[i];
         pdBaseMass[i] = a.pdBaseMass[i];
         pdNTMassUser[i] = a.pdNTMassUser[i];
      }

      return *this;
   }
};

struct ToleranceParams
{
   DoubleRange precursorTolerance;                 // raw tolerance value from param file, use to convert to true tolerace according to unit and type
   DoubleRange fragmentTolerance;                  
   int         iPrecursorToleranceUnits;           // 0=amu, 1=mmu, 2=ppm
   int         iPrecursorToleranceType;            // 0=MH+ (default), 1=precursor m/z; only valid if iMassToleranceUnits > 0
   int         iFragmentToleranceUnits;            
   int         iFragmentToleranceType;             
   int         iIsotopeSearchType;                 // different isotope search choice , for example -3, -2, -1, 0, +1, +2, +3 C13_diff

   ToleranceParams()
   {
      iPrecursorToleranceUnits = 2;
      iPrecursorToleranceType = 1;
      iFragmentToleranceUnits = 2;
      iFragmentToleranceType = 1;
      iIsotopeSearchType = 2;
      precursorTolerance = DoubleRange(3.0, 3.0);
      fragmentTolerance = DoubleRange(3.0, 3.0);
   }

   ToleranceParams& operator=(const ToleranceParams& other) {
      if (this != &other) { 
         iPrecursorToleranceUnits = other.iPrecursorToleranceUnits;
         iPrecursorToleranceType = other.iPrecursorToleranceType;
         iFragmentToleranceUnits = other.iFragmentToleranceUnits;
         iFragmentToleranceType = other.iFragmentToleranceType;
         iIsotopeSearchType = other.iIsotopeSearchType;
         precursorTolerance = other.precursorTolerance;
         fragmentTolerance = other.fragmentTolerance;
      }
      return *this;
   }

};

struct IonInfo
{
   vector<IonSeries::Ions> vSelectedIonSeries;
   int bUseWaterAmmoniaLoss;                       // ammonia, water loss
   int iTheoreticalFragmentIons;

   IonInfo() 
   {
       bUseWaterAmmoniaLoss = 0;
       iTheoreticalFragmentIons = 0;
   }

   IonInfo& operator=(const IonInfo& other)
   {
      vSelectedIonSeries = other.vSelectedIonSeries;
      bUseWaterAmmoniaLoss = other.bUseWaterAmmoniaLoss;
      iTheoreticalFragmentIons = other.iTheoreticalFragmentIons;
      return *this;
   }

};

// static user params, won't change per thread - can make global!
struct StaticParams
{
   Options              options;                // some user defined options
   DBInfo               databaseInfo;           // fasta name && info
   ToleranceParams      tolerances;
   ModificationInfo     modificationInfo;
   Enzyme               enzyme;                 // 理论酶解的酶信息, 改成只有一个，若需要成为多个，改成vector
   MassUtil             massUtility;
   IonInfo              ionInformation;
   double               dInverseBinWidth;       // this is used in BIN() many times so use inverse binWidth to do multiply vs. divide
   double               dOneMinusBinOffset;     // this is used in BIN() many times so calculate once
   int                  iXcorrProcessingOffset;

   StaticParams()
   {
      options = Options();
      databaseInfo = DBInfo();
      tolerances = ToleranceParams();
      modificationInfo = ModificationInfo();
      enzyme = Enzyme();
      massUtility = MassUtil();
      ionInformation = IonInfo();
      dInverseBinWidth = 50;
      dOneMinusBinOffset = 0.0;
      iXcorrProcessingOffset = 75;
   }

   StaticParams& operator=(const StaticParams& a)
   {
      options = a.options;
      databaseInfo = a.databaseInfo;
      tolerances = a.tolerances;
      modificationInfo = a.modificationInfo;
      enzyme = a.enzyme;
      massUtility = a.massUtility;
      ionInformation = a.ionInformation;
      dInverseBinWidth = a.dInverseBinWidth;
      dOneMinusBinOffset = a.dOneMinusBinOffset;
      iXcorrProcessingOffset = a.iXcorrProcessingOffset;
      return *this;
   }

};


class Param
{
public:
    static StaticParams                                 g_staticParams;         //static paramter, should not change in the program
    static string                                       g_sNebulionVersion;
    static MassChargeRange                              g_mcRange;
    static vector<InputFile>                            g_vInputFiles;
    static vector<NucleicAcid>                          g_vNucleicAcid;
    static vector<Oligonucleotide>                      g_vOligonucleotide;
    //static Status                                       g_cometStatus;
    static vector<VarModification>                      g_vVarModification;
    static Mutex                                        g_mcRangeMutex;
    static Mutex                                        g_vInputFileMutex;
    static Mutex                                        g_vNucleicAcidMutex;
    static Mutex                                        g_vOligonucleotideMutex;
    static Mutex                                        g_vVarModificationMutex;

    static void initMutex();

    static void destroyMutex();

    static bool loadParam(string sParamFilename);

    static void printParams(string sParamFilename, int option);

    static bool ProcessCmdLine(int argc, char* argv[]);

    static bool printUsage(string sCmd);

    static bool CheckParamsValidation();

private:

    static bool SetOptions(string arg);

};



#endif
