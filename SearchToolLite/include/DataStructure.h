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

   string toString() const
   {
      return to_string(dStart) + " -- " + to_string(dEnd);
   }

   bool inRange(double value)
   {
       if (this->dStart <= value && value <= this->dEnd)
           return true;
       return false;
   }

   static bool inRange(DoubleRange range, double value)
   {
      if (range.dStart <= value && value <= range.dEnd)
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

   IntRange& operator*(const int& a)
   {
      iStart *= a;
      iEnd *= a;
      return *this;
   }

   IntRange& operator+(const int& a)
   {
      iStart += a;
      iEnd += a;
      return *this;
   }

   string toString() const
   {
      return to_string(iStart) + " -- " + to_string(iEnd);
   }

   bool inRange(int value)
   {
       if (this->iStart <= value && value <= this->iEnd)
           return true;
       return false;
   }

   static bool inRange(IntRange range, int value)
   {
      if (range.iStart <= value && value <= range.iEnd)
         return true;
      return false;
   }
};

struct Scores
{
    double xCorr;
    double dSp;
    double dCn;
    double dExpect;
    int matchedIons;
    int totalIons;

    Scores()
    {
      xCorr = dSp = dCn = dExpect = matchedIons = totalIons = 0;
    }

    Scores(const Scores& other)
    {
      xCorr = other.xCorr;
      dSp = other.dSp;
      dCn = other.dCn;
      dExpect = other.dExpect;
      matchedIons = other.matchedIons;
      totalIons = other.totalIons;
    }

    Scores& operator=(const Scores& a)
    {
      if (this != &a)
      {
         xCorr = a.xCorr;
         dSp = a.dSp;
         dCn = a.dCn;
         dExpect = a.dExpect;
         matchedIons = a.matchedIons;
         totalIons = a.totalIons;
      }
        return *this;
    }
};

struct ExpSpectrum
{
   Spectrum spectrum;
   double   dExpPreMass;                                 //M mass
   int      iPreCharge;
   DoubleRange preBasicTolerance;
   DoubleRange preMassSearchRange;
   vector<DoubleRange> vPreMassSearchRange;      
   int      iArraySize;                                  // m/z versus intensity array
   int      iHighestIonIndex;
   double   dTotalIntensity;

   ExpSpectrum()
   {
       dExpPreMass = 0.0;
       iPreCharge = 0;
       preBasicTolerance = DoubleRange(0.0, 0.0);
       preMassSearchRange = DoubleRange(0.0, 0.0);
       iArraySize = 0;
       iHighestIonIndex = 0;
       dTotalIntensity = 0.0;
   }

   ExpSpectrum(const ExpSpectrum& other)
   {
      spectrum = other.spectrum;
      dExpPreMass = other.dExpPreMass;
      iPreCharge = other.iPreCharge;
      preBasicTolerance = other.preBasicTolerance;
      preMassSearchRange = other.preMassSearchRange;
      iArraySize = other.iArraySize;
      iHighestIonIndex = other.iHighestIonIndex;
      dTotalIntensity = other.dTotalIntensity;
      vPreMassSearchRange = other.vPreMassSearchRange;
   }

   ExpSpectrum& operator=(const ExpSpectrum& other)
   {
      if (this != &other) 
      {
         spectrum = other.spectrum;
         dExpPreMass = other.dExpPreMass;
         iPreCharge = other.iPreCharge;
         preBasicTolerance = other.preBasicTolerance;
         preMassSearchRange = other.preMassSearchRange;
         iArraySize = other.iArraySize;
         iHighestIonIndex = other.iHighestIonIndex;
         dTotalIntensity = other.dTotalIntensity;
         vPreMassSearchRange = other.vPreMassSearchRange;
      }
      return *this;
   }

};

struct IonSeries
{
   enum Ions
   {
      ion_a = 0,
      ion_b,
      ion_c,
      ion_d,
      ion_w,
      ion_x,
      ion_y,
      ion_z,
      ion_a_B,
      ion_NL,
   };
   static constexpr int iIonSeriesNum = 10;

   static double getTermMassMod(Ions ionType, bool free = false)
   {
      double dModMass = 0.0;
      double C = ElementMassUtils::Carbon_Mono;
      double H = ElementMassUtils::Hydrogen_Mono;
      double O = ElementMassUtils::Oxygen_Mono;
      double P = ElementMassUtils::Phosphorus_Mono;
      switch (ionType)
      {
      case ion_a:
      case ion_a_B:
         dModMass -= (H + 3 * O + P);
         if (!free)
            dModMass -= H;
         return dModMass;

      case ion_b:
         dModMass -= (H + 2 * O + P);
         if (!free)
            dModMass += H;
         return dModMass;

      case ion_c:
         dModMass -= 0.0;
         if (!free)
            dModMass -= H;
         return dModMass;

      case ion_d:
         dModMass += O;
         if (!free)
            dModMass += H;
         return dModMass;
      
      case ion_w:
         dModMass += (H + 3 * O + P);
         if (!free)
            dModMass += H;
         return dModMass;
      
      case ion_x:
         dModMass += (H + 2 * O + P);
         if (!free)
            dModMass -= H;
         return dModMass;

      case ion_y:
         dModMass -= 0.0;
         if (!free)
            dModMass += H;
         return dModMass;

      case ion_z:
         dModMass -= O;
         if (!free)
            dModMass -= H;
         return dModMass;

      default: //NL_ion
         return 0;
      }
   }

   static int getBYTerm(Ions ionType) // true = B-term ; false = Y-term
   {
      switch (ionType)
      {
      case ion_a:
      case ion_a_B:
      case ion_b:
      case ion_c:
      case ion_d:
         return 1;
         break;
      case ion_w:
      case ion_x:
      case ion_y:
      case ion_z:
         return 0;
         break;
      
      default: //NL_ion
         return -1;
         break;
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

struct VarModification
{
   string sSequence;
   int iMaxModificationNumber;                                                // the max modification NT in sSequence, if sSequnce size is 4, but the max number is 3, each various modification only modify 3 nt
   vector<vector<pair<int, Modification>>>   varModificationList;             // all possible of various modification, int is the index of modification
   vector<double>                            vVarModificationMassList;
   DoubleRange                               varModMassRange;

   VarModification() 
      : sSequence(""), iMaxModificationNumber(0) { }

   VarModification(const VarModification& other)
      : sSequence(other.sSequence),
         iMaxModificationNumber(other.iMaxModificationNumber),
         varModificationList(other.varModificationList),
         vVarModificationMassList(other.vVarModificationMassList),
         varModMassRange(other.varModMassRange) {}
   VarModification& operator=(const VarModification& other) 
   {
      if (this != &other) {
         sSequence = other.sSequence;
         iMaxModificationNumber = other.iMaxModificationNumber;
         varModificationList = other.varModificationList;
         vVarModificationMassList = other.vVarModificationMassList;
         varModMassRange = other.varModMassRange;
      }
      return *this;
   }

   void CalcMassRange()
   {
      varModMassRange.dStart = varModMassRange.dEnd = 0.0;
      if(vVarModificationMassList.size() == 0 || varModificationList.size() == 0)
      {
         return;
      }

      varModMassRange.dStart = *min_element(vVarModificationMassList.begin(), vVarModificationMassList.end());
      varModMassRange.dEnd = *max_element(vVarModificationMassList.begin(), vVarModificationMassList.end());
   }
};

struct NucleicAcid
{
   string sName;
   string sSequence;
   unsigned long long lFilePosition;
   int  iWhichNucleicAcid;                                     // the index in the vector of nucleic acid, useless?
   vector<pair<int, Modification>> vMods;                      // the list of Static Modification of Nucleic Acid
   Modification end5TermMod;
   Modification end3TermMod;

   NucleicAcid()
   {
      sName = "";
      sSequence = "";
      lFilePosition = 0;
      iWhichNucleicAcid = 0;
      end5TermMod = Modification();
      end3TermMod = Modification();
   }

   NucleicAcid(const NucleicAcid& other)
   {
      sName = other.sName;
      sSequence = other.sSequence;
      lFilePosition = other.lFilePosition;
      iWhichNucleicAcid = other.iWhichNucleicAcid;
      end5TermMod = other.end5TermMod;
      end3TermMod = other.end3TermMod;
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
         end5TermMod = other.end5TermMod;
         end3TermMod = other.end3TermMod;
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
   double dMass;
   string sSequence;                                                 // the basic represent of fragment, from 5' to 3', AGCTU etc.
   vector<NucleicAcidIndex> whichNucleicAcid;
   bool bDecoy;                                                      // is the oligonucleotide is decoys
   vector<pair<int, Modification>> vMods;                            // the list of **static** Modifications on each position
   int iVarModificationIndex;
   Modification end5TermMod;
   Modification end3TermMod;

   Oligonucleotide()
   {
      dMass = 0.0;
      whichNucleicAcid = vector<NucleicAcidIndex>();
      sSequence = "";
      iVarModificationIndex = 0;
      end5TermMod = Modification();
      end3TermMod = Modification();
      bDecoy = false;
   }

   Oligonucleotide(const Oligonucleotide& other)
   {
      dMass = other.dMass;
      whichNucleicAcid = other.whichNucleicAcid;
      sSequence = other.sSequence;
      iVarModificationIndex = other.iVarModificationIndex;
      vMods = other.vMods;
      end5TermMod = other.end5TermMod;
      end3TermMod = other.end3TermMod;
      bDecoy = other.bDecoy;
   }

   Oligonucleotide& operator=(const Oligonucleotide& other)
   {
      if (this != &other) { 
            dMass = other.dMass;
            whichNucleicAcid = other.whichNucleicAcid;
            sSequence = other.sSequence;
            iVarModificationIndex = other.iVarModificationIndex;
            vMods = other.vMods;
            end5TermMod = other.end5TermMod;
            end3TermMod = other.end3TermMod;
            bDecoy = other.bDecoy;
      }
      return *this;  
   }

};

struct OligonucleotideFragment
{
   double dMass;
   IonSeries::Ions ionType;
   string sSequence;                                                 // the basic represent of fragment, from 5' to 3', AGCTU etc.
   int iCharge;
   vector<pair<int, Modification>> vMods;                            // the list of **all** Modifications on each position
   Modification end5TermMod;
   Modification end3TermMod;

   OligonucleotideFragment()
   {
       dMass = 0.0;
       ionType = IonSeries::ion_a;
       iCharge = 0;
       sSequence = "";
       end5TermMod = Modification();
       end3TermMod = Modification();
   }

   OligonucleotideFragment(const OligonucleotideFragment& other)
   {
       dMass = other.dMass;
       ionType = other.ionType;
       iCharge = other.iCharge;
       sSequence = other.sSequence;
       vMods = other.vMods;
       end5TermMod = other.end5TermMod;
       end3TermMod = other.end3TermMod;
   }

   OligonucleotideFragment& operator=(const OligonucleotideFragment& other) 
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
   struct QueryEntry
   {
      Oligonucleotide&      oligo;
      int                         iWhichVarModCombination;
      Scores                      scores;

      QueryEntry(Oligonucleotide& oligo): oligo(oligo)
      {
         iWhichVarModCombination = 0;
         scores = Scores();
      }

      QueryEntry(const QueryEntry& other): oligo(other.oligo)
      {
         iWhichVarModCombination = other.iWhichVarModCombination;
         scores = other.scores;
      }

      QueryEntry& operator=(const QueryEntry other)
      {
          if (this != &other)
          {
              this->oligo = other.oligo;
              this->iWhichVarModCombination = other.iWhichVarModCombination;
              this->scores = other.scores;
          }
          return *this;
      }

   };

   vector<QueryEntry> vResults;
   ExpSpectrum*       expSpectrum;                              // experiment spectrum 

   Mutex accessMutex;

   Query()
   {
       expSpectrum = nullptr;
       Threading::CreateMutex(&accessMutex);
   }

   ~Query()
   {
       Threading::DestroyMutex(accessMutex);
   }

   void pushBackResult(Oligonucleotide& oligo, int iWhichVarModCombination)
   {
      QueryEntry entry(oligo);
      entry.iWhichVarModCombination = iWhichVarModCombination;
      vResults.push_back(entry);
   }

};

#endif // _COMETDATA_H_
