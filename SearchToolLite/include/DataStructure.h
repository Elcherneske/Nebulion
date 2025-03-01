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


#ifndef _DATASTRUCTURE_H_
#define _DATASTRUCTURE_H_

#include "Common.h"
#include "Threading.h"

#include <vector>
#include <string>

struct DoubleRange
{
   double dStart;
   double dEnd;

   DoubleRange()
   {
      dStart = 0.0;
      dEnd = 0.0;
   }

   DoubleRange(const DoubleRange& a)
   {
      dStart = a.dStart;
      dEnd = a.dEnd;
   }

   DoubleRange(double dStart_in, double dEnd_in)
   {
      dStart = dStart_in;
      dEnd = dEnd_in;
   }

   DoubleRange& operator=(const DoubleRange& a)
   {
      dStart = a.dStart;
      dEnd = a.dEnd;
      return *this;
   }

   DoubleRange& operator*=(const double& a)
   {
      dStart *= a;
      dEnd *= a;
      return *this;
   }

   DoubleRange& operator+=(const double& a)
   {
      dStart += a;
      dEnd += a;
      return *this;
   }

   bool inRange(double value) const
   {
       if (this->dStart <= value && value <= this->dEnd)
           return true;
       return false;
   }
};

struct IntRange
{
   int iStart;
   int iEnd;

   IntRange()
   {
      iStart = 0;
      iEnd = 0;
   }

   IntRange(const IntRange& a)
   {
      iStart = a.iStart;
      iEnd = a.iEnd;
   }

   IntRange(int iStart_in, int iEnd_in)
   {
      iStart = iStart_in;
      iEnd = iEnd_in;
   }

   IntRange& operator=(const IntRange& a)
   {
      iStart = a.iStart;
      iEnd = a.iEnd;
      return *this;
   }

   IntRange& operator*=(const int& a)
   {
      iStart *= a;
      iEnd *= a;
      return *this;
   }

   IntRange& operator+=(const int& a)
   {
      iStart += a;
      iEnd += a;
      return *this;
   }

   bool inRange(int value) const
   {
       if (this->iStart <= value && value <= this->iEnd)
           return true;
       return false;
   }
};

struct RoughScore
{
    double xCorr;
    double dSp;
    double dCn;
    int matchedIons;
    int totalIons;

    RoughScore() {
      xCorr = dSp = dCn = matchedIons = totalIons = 0;
    }

    RoughScore(const RoughScore& other) {
      xCorr = other.xCorr;
      dSp = other.dSp;
      dCn = other.dCn;
      matchedIons = other.matchedIons;
      totalIons = other.totalIons;
    }

    RoughScore& operator=(const RoughScore& other) {
      if (this != &other) {
         xCorr = other.xCorr;
         dSp = other.dSp;
         dCn = other.dCn;
         matchedIons = other.matchedIons;
         totalIons = other.totalIons;
      }
        return *this;
    }
};

struct ExpSpectrum
{
   struct Peak
   {
      double dMZ;
      double dIntensity;

      Peak(double dMZ, double dIntensity)
      {
         this->dMZ = dMZ;
         this->dIntensity = dIntensity;
      }

      Peak& operator=(const Peak& other)
      {
         if (this != &other)
         {
            dMZ = other.dMZ;
            dIntensity = other.dIntensity;
         }
         return *this;
      }
   };
   std::vector<Peak> spectrum;
   double   dExpMZ; //experiment m/z
   int      iPreCharge;
   DoubleRange preMassSearchRange;
   vector<DoubleRange> vPreMassMatchRange;      
   int      iArraySize;// m/z versus intensity array
   int      iHighestIonIndex;
   double   dTotalIntensity;

   ExpSpectrum()
   {
      spectrum = vector<Peak>();
      dExpMZ = 0.0;
      iPreCharge = 0;
      preMassSearchRange = DoubleRange(0.0, 0.0);
      vPreMassMatchRange = vector<DoubleRange>();
      iArraySize = 0;
      iHighestIonIndex = 0;
      dTotalIntensity = 0.0; 
   }

   ExpSpectrum(const ExpSpectrum& other)
   {
      spectrum = other.spectrum;
      dExpMZ = other.dExpMZ;
      iPreCharge = other.iPreCharge;
      preMassSearchRange = other.preMassSearchRange;
      iArraySize = other.iArraySize;
      iHighestIonIndex = other.iHighestIonIndex;
      dTotalIntensity = other.dTotalIntensity;
      vPreMassMatchRange = other.vPreMassMatchRange;
   }

   ExpSpectrum& operator=(const ExpSpectrum& other)
   {
      if (this != &other) 
      {
         spectrum = other.spectrum;
         dExpMZ = other.dExpMZ;
         iPreCharge = other.iPreCharge;
         preMassSearchRange = other.preMassSearchRange;
         iArraySize = other.iArraySize;
         iHighestIonIndex = other.iHighestIonIndex;
         dTotalIntensity = other.dTotalIntensity;
         vPreMassMatchRange = other.vPreMassMatchRange;
      }
      return *this;
   }

};

struct IonSpecies
{
   enum IonID
   {
      ion_a = 0,
      ion_b,
      ion_c,
      ion_d,
      ion_w,
      ion_x,
      ion_y,
      ion_z,
      ion_NL,
   };
   static constexpr int iIonSeriesNum = 9;

   IonID ionType;
   bool bFreeRadical;
   bool bDropBase;

   IonSpecies()
   {
      ionType = ion_a;
      bFreeRadical = false;
      bDropBase = false;
   }

   IonSpecies(IonID ionType, bool bFreeRadical, bool bDropBase)
   {
      this->ionType = ionType;
      this->bFreeRadical = bFreeRadical;
      this->bDropBase = bDropBase;
   }

   IonSpecies(const IonSpecies& other)
   {
      ionType = other.ionType;
      bFreeRadical = other.bFreeRadical;
      bDropBase = other.bDropBase;
   }

   IonSpecies& operator=(const IonSpecies& other)
   {
      if (this != &other)
      {
         ionType = other.ionType;
         bFreeRadical = other.bFreeRadical;
         bDropBase = other.bDropBase;
      }
      return *this;
   }

   int BYIonType() // 1 = B-term ; 0 = Y-term
   {
      switch (this->ionType)
      {
         case ion_a:
         case ion_b:
         case ion_c:
         case ion_d:
            return 1;
         case ion_w:
         case ion_x:
         case ion_y:
         case ion_z:
            return 0;
         default: //NL_ion
            return -1;
      }
   }
};

struct Modification
{
   double dModMass;
   string sModName;
   int iModID;
   
   Modification()
   {
      dModMass = 0.0;
      sModName = "";
      iModID = 0;
   }

   Modification(const Modification& other)
   {
      dModMass = other.dModMass;
      sModName = other.sModName;
      iModID = other.iModID;
   }

   Modification& operator=(const Modification& other) 
   {
      if (this != &other) {
         dModMass = other.dModMass;
         sModName = other.sModName;
         iModID = other.iModID;
      }
      return *this;
   }

   string toString() const
   {
      return "[" + to_string(iModID) + ";" + sModName + ";" + to_string(dModMass) + "]";
   }

};

struct VarModFormula
{
   string sFormula;
   vector<vector<Modification>>   vvVarModList; // all possible of various modification, int is the index of modification
   vector<double>                 vVarModMassList;
   DoubleRange                    varModMassRange;

   VarModFormula()
   {
      sFormula = "";
      vvVarModList = vector<vector<Modification>>();
      vVarModMassList = vector<double>();
      varModMassRange = DoubleRange();
   }

   VarModFormula(const VarModFormula& other)
   {
      sFormula = other.sFormula;
      vvVarModList = other.vvVarModList;
      vVarModMassList = other.vVarModMassList;
      varModMassRange = other.varModMassRange;
   }

   VarModFormula& operator=(const VarModFormula& other) 
   {
      if (this != &other) {
         sFormula = other.sFormula;
         vvVarModList = other.vvVarModList;
         vVarModMassList = other.vVarModMassList;
         varModMassRange = other.varModMassRange;
      }
      return *this;
   }

   void CalcMassRange()
   {
      varModMassRange.dStart = varModMassRange.dEnd = 0.0;
      if(vVarModMassList.size() == 0 || vvVarModList.size() == 0)
      {
         return;
      }
      varModMassRange.dStart = *min_element(vVarModMassList.begin(), vVarModMassList.end());
      varModMassRange.dEnd = *max_element(vVarModMassList.begin(), vVarModMassList.end());
   }
};

struct NucleicAcid
{
   string sName;
   string sSequence;
   unsigned long long lFilePosition;
   int  iWhichNucleicAcid;                                     // the index in the vector of nucleic acid, useless?
   vector<pair<int, Modification>> vMods;                      // the list of Static Modification of Nucleic Acid

   NucleicAcid()
   {
      sName = "";
      sSequence = "";
      lFilePosition = 0;
      iWhichNucleicAcid = 0;
   }

   NucleicAcid(const NucleicAcid& other)
   {
      sName = other.sName;
      sSequence = other.sSequence;
      lFilePosition = other.lFilePosition;
      iWhichNucleicAcid = other.iWhichNucleicAcid;
      vMods = other.vMods;
   }

   NucleicAcid& operator=(const NucleicAcid& other) 
   {
      if (this != &other) { 
         sName = other.sName;
         sSequence = other.sSequence;
         lFilePosition = other.lFilePosition;
         iWhichNucleicAcid = other.iWhichNucleicAcid;
         vMods = other.vMods;
      }
      return *this;
   }

};

struct Oligonucleotide
{
   struct NucleicAcidIndex
   {
      int iWhichNucleicAcid;
      int iNABeginPos;
      int iNAEndPos;

      NucleicAcidIndex()
      {
         iWhichNucleicAcid = iNABeginPos = iNAEndPos = 0;
      }

      NucleicAcidIndex(const NucleicAcidIndex& other)
      {
         iWhichNucleicAcid = other.iWhichNucleicAcid;
         iNABeginPos = other.iNABeginPos;
         iNAEndPos = other.iNAEndPos;
      }

      NucleicAcidIndex& operator=(const NucleicAcidIndex& other)
      {
         if (this != &other) { 
               iWhichNucleicAcid = other.iWhichNucleicAcid;
               iNABeginPos = other.iNABeginPos;
               iNAEndPos = other.iNAEndPos;
         }
         return *this;  
      }
   };
   
   string sSequence;                                                 // the basic represent of fragment, from 5' to 3', AGCTU etc.
   vector<NucleicAcidIndex> vWhichNucleicAcid;
   bool bDecoy;                                                      // is the oligonucleotide is decoys
   vector<pair<int, Modification>> vMods;                      // the list of **static** Modifications on each position
   vector<int> vVarModFormulaIndices;
   map<int, vector<vector<int>>> mVarModPosition; //{index of vVarModFormulaIndices, multiple position indices}
   Modification end5TermMod;
   Modification end3TermMod;

   Oligonucleotide()
   {
      vWhichNucleicAcid = vector<NucleicAcidIndex>();
      sSequence = "";
      vMods = vector<pair<int, Modification>>();
      vVarModFormulaIndices = vector<int>();
      mVarModPosition = map<int, vector<vector<int>>>();
      end5TermMod = Modification();
      end3TermMod = Modification();
      bDecoy = false;
   }

   Oligonucleotide(const Oligonucleotide& other)
   {
      vWhichNucleicAcid = other.vWhichNucleicAcid;
      sSequence = other.sSequence;
      vVarModFormulaIndices = other.vVarModFormulaIndices;
      mVarModPosition = other.mVarModPosition;
      vMods = other.vMods;
      end5TermMod = other.end5TermMod;
      end3TermMod = other.end3TermMod;
      bDecoy = other.bDecoy;
   }

   Oligonucleotide& operator=(const Oligonucleotide& other)
   {
      if (this != &other) { 
         vWhichNucleicAcid = other.vWhichNucleicAcid;
         sSequence = other.sSequence;
         vVarModFormulaIndices = other.vVarModFormulaIndices;
         mVarModPosition = other.mVarModPosition;
         vMods = other.vMods;
         end5TermMod = other.end5TermMod;
         end3TermMod = other.end3TermMod;
         bDecoy = other.bDecoy;
      }
      return *this;  
   }

};

struct OligonucleotideID
   {
      int iOligoIndex;
      int iVarModFormulaIndex;
      int iPositionIndex;
      int iVarModIndex;

      OligonucleotideID()
      {
         iOligoIndex = iVarModFormulaIndex = iPositionIndex = iVarModIndex = 0;
      }

      OligonucleotideID(int iOligoIndex, int iPositionIndex, int iVarModFormulaIndex, int iVarModIndex)
      {
         this->iOligoIndex = iOligoIndex;
         this->iPositionIndex = iPositionIndex;
         this->iVarModFormulaIndex = iVarModFormulaIndex;
         this->iVarModIndex = iVarModIndex;
      }

      OligonucleotideID(const OligonucleotideID& other)
      {
         iOligoIndex = other.iOligoIndex;
         iVarModFormulaIndex = other.iVarModFormulaIndex;
         iPositionIndex = other.iPositionIndex;
         iVarModIndex = other.iVarModIndex;
      }

      OligonucleotideID& operator=(const OligonucleotideID& other)
      {
         if (this != &other)
         {
            iOligoIndex = other.iOligoIndex;
            iVarModFormulaIndex = other.iVarModFormulaIndex;
            iPositionIndex = other.iPositionIndex;
            iVarModIndex = other.iVarModIndex;
         }
         return *this;
      }

      bool operator==(const OligonucleotideID& other) const
      {
         return 
            iOligoIndex == other.iOligoIndex && 
            iVarModFormulaIndex == other.iVarModFormulaIndex && 
            iPositionIndex == other.iPositionIndex && 
            iVarModIndex == other.iVarModIndex;
      }

      bool operator<(const OligonucleotideID& other) const
      {
         if (iOligoIndex != other.iOligoIndex)
            return iOligoIndex < other.iOligoIndex;
         if (iPositionIndex != other.iPositionIndex)
            return iPositionIndex < other.iPositionIndex;
         if (iVarModFormulaIndex != other.iVarModFormulaIndex)
            return iVarModFormulaIndex < other.iVarModFormulaIndex;
         return iVarModIndex < other.iVarModIndex;
      }
   };

struct OligoFragment
{
   double dMass;
   IonSpecies ionType;
   std::string sSequence;                                                 // the basic represent of fragment, from 5' to 3', AGCTU etc.
   int iCharge;
   vector<pair<int, Modification>> vMods;                            // the list of **all** Modifications on each position
   Modification end5TermMod;
   Modification end3TermMod;

   OligoFragment()
   {
       dMass = 0.0;
       ionType = IonSpecies();
       iCharge = 0;
       sSequence = "";
       end5TermMod = Modification();
       end3TermMod = Modification();
   }

   OligoFragment(const OligoFragment& other)
   {
       dMass = other.dMass;
       ionType = other.ionType;
       iCharge = other.iCharge;
       sSequence = other.sSequence;
       vMods = other.vMods;
       end5TermMod = other.end5TermMod;
       end3TermMod = other.end3TermMod;
   }

   OligoFragment& operator=(const OligoFragment& other) 
   {
      if (this != &other) {
         dMass = other.dMass;
         ionType = other.ionType;
         iCharge = other.iCharge;
         sSequence = other.sSequence;
         vMods = other.vMods;
         end5TermMod = other.end5TermMod;
         end3TermMod = other.end3TermMod;
      }
      return *this;
   }
};

struct Enzyme
{
   int  iAllowedMissedCleavage;
   int  iSearchEnzymeOffSet;
   int  iEnzymeSelected;
   string sEnzymeName;
   string sEnzymeBreakNT;
   string sEnzymeNoBreakNT;

   Enzyme(): iAllowedMissedCleavage(0), iEnzymeSelected(0), iSearchEnzymeOffSet(0), sEnzymeName(""), sEnzymeBreakNT(""), sEnzymeNoBreakNT("") { }

   Enzyme(const Enzyme& other)
      : iAllowedMissedCleavage(other.iAllowedMissedCleavage), iEnzymeSelected(other.iEnzymeSelected), iSearchEnzymeOffSet(other.iSearchEnzymeOffSet), 
      sEnzymeName(other.sEnzymeName), sEnzymeBreakNT(other.sEnzymeBreakNT), sEnzymeNoBreakNT(other.sEnzymeNoBreakNT) { }

   Enzyme& operator=(const Enzyme& other) 
   {
      if (this != &other) {  
         iAllowedMissedCleavage = other.iAllowedMissedCleavage;
         iSearchEnzymeOffSet = other.iSearchEnzymeOffSet;
         iEnzymeSelected = other.iEnzymeSelected;
         sEnzymeName = other.sEnzymeName;
         sEnzymeBreakNT = other.sEnzymeBreakNT;
         sEnzymeNoBreakNT = other.sEnzymeNoBreakNT;
      }
      return *this;
   }

   bool operator==(const Enzyme& other) const 
   {
      if (iAllowedMissedCleavage != other.iAllowedMissedCleavage) return false;
      if (iSearchEnzymeOffSet != other.iSearchEnzymeOffSet) return false;
      if (iEnzymeSelected != other.iEnzymeSelected) return false;
      if (sEnzymeName != other.sEnzymeName) return false;
      if (sEnzymeBreakNT != other.sEnzymeBreakNT) return false;
      if (sEnzymeNoBreakNT != other.sEnzymeNoBreakNT) return false;
      return true;
   }

};

struct InputFile
{
   enum AnalysisType
   {
      AnalysisType_Unknown = 0,
      AnalysisType_DTA,
      AnalysisType_SpecificScan,
      AnalysisType_SpecificScanRange,
      AnalysisType_EntireFile
   };
   
   enum InputType
   {
      InputType_UNKNOWN = -1,
      InputType_MS2 = 0,           // ms2, cms2, bms2, etc.
      InputType_MZXML,
      InputType_MZML,
      InputType_RAW,
      InputType_MGF
   };

   int  iInputType;           // reference to enum inputType
   int  iAnalysisType;        // reference to enum AnalysisType
   int  iFirstScan;           // for scan range specified with file on command line; otherwise replicates g_staticParams.options.scanRange.iStart
   int  iLastScan;
   string sFileName;          //filename
   string sBaseName;          //the database name

   InputFile()
   {
      iInputType = 0;
      iAnalysisType = AnalysisType_Unknown;
      iFirstScan = 0;
      iLastScan = 0;
      sFileName = "";
      sBaseName = "";
   }

   InputFile(const InputFile& inputObj)
   {
      iInputType = inputObj.iInputType;
      iAnalysisType = inputObj.iAnalysisType;
      iFirstScan = inputObj.iFirstScan;
      iLastScan = inputObj.iLastScan;
      sFileName = inputObj.sFileName;
      sBaseName = inputObj.sBaseName;
   }

   InputFile& operator = (const InputFile &inputObj)
   {
      iInputType = inputObj.iInputType;
      iAnalysisType = inputObj.iAnalysisType;
      iFirstScan = inputObj.iFirstScan;
      iLastScan = inputObj.iLastScan;
      sBaseName = inputObj.sBaseName;
      sFileName = inputObj.sFileName;
      return *this;
   }
};

struct Query
{
   struct OligoEntry
   {
      OligonucleotideID          oligoID;
      RoughScore                  scores;

      OligoEntry(int iOligoIndex, int iPositionIndex, int iVarModFormulaIndex, int iVarModIndex)
      {
         oligoID = OligonucleotideID(iOligoIndex, iPositionIndex, iVarModFormulaIndex, iVarModIndex);
         scores = RoughScore();
      }
      
      OligoEntry(const OligonucleotideID& oligoID)
      {
         this->oligoID = OligonucleotideID(oligoID);
         scores = RoughScore();
      }

      OligoEntry(const OligoEntry& other)
      {
         oligoID = other.oligoID;
         scores = other.scores;
      }
   };

   vector<OligoEntry> vOligoEntries;
   ExpSpectrum&       expSpectrum;                              // experiment spectrum 
   Query(ExpSpectrum& expSpectrum): expSpectrum(expSpectrum)
   {}

   ~Query()
   {}

};

#endif
