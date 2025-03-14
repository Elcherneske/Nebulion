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

#include "OligoProcessUtils.h"
#include "search.h"
#include "Status.h"
#include "MassUtils.h"
#include "ModificationUtils.h"
#include "WriteOutUtils.h"

#include <map>
#include <set>
#include <sstream>
#include <fstream>
#include <random>

OligoProcessUtils::OligoProcessUtils()
{
}

OligoProcessUtils::~OligoProcessUtils()
{
}

vector<OligoFragment> OligoProcessUtils::GenerateFragments(const Oligonucleotide& oligo, int iCharge, const OligonucleotideID& oligoID)
{
   vector<OligoFragment> vFragmentResults;
   if (Param::getInstance().options.bUseNLIons)
      vFragmentResults = OligoProcessUtils::GenerateNLFragments(oligo, iCharge, oligoID);

   int iOligoLen = oligo.sSequence.size();
   int iPosReverse;                       // points to residue in reverse order

   if (iOligoLen <= 2)
      return vFragmentResults;
   
   for (int iIndex = 1; iIndex < oligo.sSequence.size() - 1; ++iIndex) // fragment length is larger than 2
   {
      iPosReverse = oligo.sSequence.size() - iIndex - 1;
      int IonLen = iIndex + 1;

      for (int iIonCharge = 1; iIonCharge <= abs(iCharge); iIonCharge++)
      {
         OligoFragment bTermIon;
         OligoFragment yTermIon;

         bTermIon.sSequence = oligo.sSequence.substr(0, IonLen);
         bTermIon.end5TermMod = oligo.end5TermMod;
         if (Param::getInstance().options.bPositiveCharge)
            bTermIon.iCharge = iIonCharge;
         else
            bTermIon.iCharge = -1 * iIonCharge;
         // add static modification
         for (const auto& mod : oligo.vMods) {
            int iModPos = mod.first;
            if (iModPos >= 0 && iModPos < IonLen)
               bTermIon.vMods.push_back({iModPos, mod.second});
         }
         // add variable modification
         vector<Modification> vVarMods = Param::getInstance().vVarModFormulas[oligoID.iVarModFormulaIndex].vvVarModList[oligoID.iVarModIndex];
         vector<int> vVarModIndex = oligo.mVarModPosition.at(oligoID.iVarModFormulaIndex)[oligoID.iPositionIndex];
         for (size_t iIndex = 0; iIndex < vVarModIndex.size(); iIndex++) {
            if (vVarModIndex[iIndex] >= 0 && vVarModIndex[iIndex] < IonLen)
               bTermIon.vMods.push_back({vVarModIndex[iIndex], vVarMods[iIndex]});
         }
         // consider different ion types -> 3' modification
         if (Param::getInstance().options.bUseAIons){
            Modification cutMod;
            IonSpecies ionSpecies;
            ionSpecies.ionType = IonSpecies::ion_a;
            ionSpecies.bFreeRadical = Param::getInstance().options.bFreeRadical;
            ionSpecies.bDropBase = false;
            cutMod.sModName = "a ion modification";
            cutMod.dModMass = MassUtils::getInstance().CalCutModMass(ionSpecies, bTermIon.sSequence[bTermIon.sSequence.size() - 1]);
            bTermIon.end3TermMod = cutMod;
            bTermIon.dMass = MassUtils::getInstance().CalFragmentMass(bTermIon, true, true);
            vFragmentResults.push_back(bTermIon);
         }

         if (Param::getInstance().options.bUseABIons){
            Modification cutMod;
            IonSpecies ionSpecies;
            ionSpecies.ionType = IonSpecies::ion_a;
            ionSpecies.bFreeRadical = Param::getInstance().options.bFreeRadical;
            ionSpecies.bDropBase = true;
            cutMod.sModName = "a ion modification";
            cutMod.dModMass = MassUtils::getInstance().CalCutModMass(ionSpecies, bTermIon.sSequence[bTermIon.sSequence.size() - 1]);
            bTermIon.end3TermMod = cutMod;
            bTermIon.dMass = MassUtils::getInstance().CalFragmentMass(bTermIon, true, true);
            vFragmentResults.push_back(bTermIon);
         }

         if (Param::getInstance().options.bUseBIons){
            Modification cutMod;
            IonSpecies ionSpecies;
            ionSpecies.ionType = IonSpecies::ion_b;
            ionSpecies.bFreeRadical = Param::getInstance().options.bFreeRadical;
            ionSpecies.bDropBase = false;
            cutMod.sModName = "b ion modification";
            cutMod.dModMass = MassUtils::getInstance().CalCutModMass(ionSpecies, bTermIon.sSequence[bTermIon.sSequence.size() - 1]);
            bTermIon.end3TermMod = cutMod;
            bTermIon.dMass = MassUtils::getInstance().CalFragmentMass(bTermIon, true, true);
            vFragmentResults.push_back(bTermIon);
         }

         if (Param::getInstance().options.bUseCIons){
            Modification cutMod;
            IonSpecies ionSpecies;
            ionSpecies.ionType = IonSpecies::ion_c;
            ionSpecies.bFreeRadical = Param::getInstance().options.bFreeRadical;
            ionSpecies.bDropBase = false;
            cutMod.sModName = "c ion modification";
            cutMod.dModMass = MassUtils::getInstance().CalCutModMass(ionSpecies, bTermIon.sSequence[bTermIon.sSequence.size() - 1]);
            bTermIon.end3TermMod = cutMod;
            bTermIon.dMass = MassUtils::getInstance().CalFragmentMass(bTermIon, true, true);
            vFragmentResults.push_back(bTermIon);
         }

         if (Param::getInstance().options.bUseDIons){
            Modification cutMod;
            IonSpecies ionSpecies;
            ionSpecies.ionType = IonSpecies::ion_d;
            ionSpecies.bFreeRadical = Param::getInstance().options.bFreeRadical;
            ionSpecies.bDropBase = false;
            cutMod.sModName = "d ion modification";
            cutMod.dModMass = MassUtils::getInstance().CalCutModMass(ionSpecies, bTermIon.sSequence[bTermIon.sSequence.size() - 1]);
            bTermIon.end3TermMod = cutMod;
            bTermIon.dMass = MassUtils::getInstance().CalFragmentMass(bTermIon, true, true);
            vFragmentResults.push_back(bTermIon);
         }

         yTermIon.sSequence = oligo.sSequence.substr(iPosReverse, IonLen);
         yTermIon.end3TermMod = oligo.end3TermMod;
         if (Param::getInstance().options.bPositiveCharge)
            yTermIon.iCharge = iIonCharge;
         else
            yTermIon.iCharge = -1 * iIonCharge;
         // add static modification
         for (const auto& mod : oligo.vMods) {
            int iModPos = mod.first;
            if (iModPos >= iPosReverse && iModPos < oligo.sSequence.size())
               yTermIon.vMods.push_back({iModPos - iPosReverse, mod.second});
         }
         // add variable modification
         vVarMods = Param::getInstance().vVarModFormulas[oligoID.iVarModFormulaIndex].vvVarModList[oligoID.iVarModIndex];
         vVarModIndex = oligo.mVarModPosition.at(oligoID.iVarModFormulaIndex)[oligoID.iPositionIndex];
         for (size_t iIndex = 0; iIndex < vVarModIndex.size(); iIndex++) {
            if (vVarModIndex[iIndex] >= iPosReverse && vVarModIndex[iIndex] < oligo.sSequence.size())
               yTermIon.vMods.push_back({vVarModIndex[iIndex] - iPosReverse, vVarMods[iIndex]});
         }

         if (Param::getInstance().options.bUseWIons){
            Modification cutMod;
            IonSpecies ionSpecies;
            ionSpecies.ionType = IonSpecies::ion_w;
            ionSpecies.bFreeRadical = Param::getInstance().options.bFreeRadical;
            ionSpecies.bDropBase = false;
            cutMod.sModName = "w ion modification";
            cutMod.dModMass = MassUtils::getInstance().CalCutModMass(ionSpecies, yTermIon.sSequence[yTermIon.sSequence.size() - 1]);
            yTermIon.end5TermMod = cutMod;
            yTermIon.dMass = MassUtils::getInstance().CalFragmentMass(yTermIon, true, true);
            vFragmentResults.push_back(yTermIon);
         }  

         if (Param::getInstance().options.bUseXIons){
            Modification cutMod;
            IonSpecies ionSpecies;
            ionSpecies.ionType = IonSpecies::ion_x;
            ionSpecies.bFreeRadical = Param::getInstance().options.bFreeRadical;
            ionSpecies.bDropBase = false;
            cutMod.sModName = "x ion modification";
            cutMod.dModMass = MassUtils::getInstance().CalCutModMass(ionSpecies, yTermIon.sSequence[yTermIon.sSequence.size() - 1]);
            yTermIon.end5TermMod = cutMod;
            yTermIon.dMass = MassUtils::getInstance().CalFragmentMass(yTermIon, true, true);
            vFragmentResults.push_back(yTermIon);
         }
         
         if (Param::getInstance().options.bUseYIons){
            Modification cutMod;
            IonSpecies ionSpecies;
            ionSpecies.ionType = IonSpecies::ion_y;
            ionSpecies.bFreeRadical = Param::getInstance().options.bFreeRadical;
            ionSpecies.bDropBase = false;
            cutMod.sModName = "y ion modification";
            cutMod.dModMass = MassUtils::getInstance().CalCutModMass(ionSpecies, yTermIon.sSequence[yTermIon.sSequence.size() - 1]);
            yTermIon.end5TermMod = cutMod;
            yTermIon.dMass = MassUtils::getInstance().CalFragmentMass(yTermIon, true, true);
            vFragmentResults.push_back(yTermIon);
         }

         if (Param::getInstance().options.bUseZIons){
            Modification cutMod;
            IonSpecies ionSpecies;
            ionSpecies.ionType = IonSpecies::ion_z;
            ionSpecies.bFreeRadical = Param::getInstance().options.bFreeRadical;
            ionSpecies.bDropBase = false;
            cutMod.sModName = "z ion modification";
            cutMod.dModMass = MassUtils::getInstance().CalCutModMass(ionSpecies, yTermIon.sSequence[yTermIon.sSequence.size() - 1]);
            yTermIon.end5TermMod = cutMod;
            yTermIon.dMass = MassUtils::getInstance().CalFragmentMass(yTermIon, true, true);
            vFragmentResults.push_back(yTermIon);
         }
      }
   }

   return vFragmentResults;
}

vector<OligoFragment> OligoProcessUtils::GenerateNLFragments(const Oligonucleotide& oligo, int iCharge, const OligonucleotideID& oligoID)
{
    vector<OligoFragment> vFragmentResults;

    double H = ElementMassUtils::Hydrogen_Mono;
    double P = ElementMassUtils::Phosphorus_Mono;
    double O = ElementMassUtils::Oxygen_Mono;
    double dOligoMass = MassUtils::getInstance().CalOligoMass(oligo, oligoID, true, true, true);

    // loss P or H2O
   OligoFragment lossHPO3Fragment;
   lossHPO3Fragment.end3TermMod = oligo.end3TermMod;
   lossHPO3Fragment.end5TermMod = oligo.end5TermMod;
   lossHPO3Fragment.ionType = IonSpecies(IonSpecies::ion_NL, false, false);
   lossHPO3Fragment.vMods = oligo.vMods;
   lossHPO3Fragment.sSequence = oligo.sSequence;
   lossHPO3Fragment.dMass = dOligoMass - (H + P + 3 * O);
   lossHPO3Fragment.iCharge = iCharge;
   vFragmentResults.push_back(lossHPO3Fragment);

   OligoFragment lossH3PO4Fragment;
   lossH3PO4Fragment.end3TermMod = oligo.end3TermMod;
   lossH3PO4Fragment.end5TermMod = oligo.end5TermMod;
   lossH3PO4Fragment.ionType = IonSpecies(IonSpecies::ion_NL, false, false);
   lossH3PO4Fragment.vMods = oligo.vMods;
   lossH3PO4Fragment.sSequence = oligo.sSequence;
   lossH3PO4Fragment.dMass = dOligoMass - (2 * H + O);
   lossH3PO4Fragment.iCharge = iCharge;
   vFragmentResults.push_back(lossH3PO4Fragment);

   OligoFragment lossH2OFragment;
   lossH2OFragment.end3TermMod = oligo.end3TermMod;
   lossH2OFragment.end5TermMod = oligo.end5TermMod;
   lossH2OFragment.ionType = IonSpecies(IonSpecies::ion_NL, false, false);
   lossH2OFragment.vMods = oligo.vMods;
   lossH2OFragment.sSequence = oligo.sSequence;
   lossH2OFragment.dMass = dOligoMass - (3 * H + P + 4 * O);
   lossH2OFragment.iCharge = iCharge;
   vFragmentResults.push_back(lossH2OFragment);

    // base fragment and M-Base fragment
    vector<char> vBase = {'A', 'T', 'C', 'G', 'U'};
    for (char base : vBase) {
        if (oligo.sSequence.find(base) != string::npos) {
            // M-base fragment
            OligoFragment dropBaseFragment;
            dropBaseFragment.end3TermMod = oligo.end3TermMod;
            dropBaseFragment.end5TermMod = oligo.end5TermMod;
            dropBaseFragment.ionType = IonSpecies(IonSpecies::ion_NL, false, true);
            dropBaseFragment.vMods = oligo.vMods;
            dropBaseFragment.sSequence = oligo.sSequence;
            if (Param::getInstance().options.bPositiveCharge)
                dropBaseFragment.iCharge = abs(iCharge) - 1;
            else
                dropBaseFragment.iCharge = -1 * abs(iCharge) + 1;
            if (Param::getInstance().options.bFreeRadical)
                dropBaseFragment.dMass = dOligoMass - (MassUtils::getInstance().GetBaseMass(base));
            else
                dropBaseFragment.dMass = dOligoMass - (MassUtils::getInstance().GetBaseMass(base) + ElementMassUtils::Hydrogen_Mono);
            vFragmentResults.push_back(dropBaseFragment);

            // base fragment
            OligoFragment baseFragment;
            baseFragment.ionType = IonSpecies(IonSpecies::ion_NL, false, false);
            baseFragment.sSequence = "" + base;
            if (Param::getInstance().options.bFreeRadical)
                baseFragment.dMass = MassUtils::getInstance().GetBaseMass(base);
            else
                baseFragment.dMass = MassUtils::getInstance().GetBaseMass(base) + ElementMassUtils::Hydrogen_Mono;
            if (Param::getInstance().options.bPositiveCharge)
                baseFragment.iCharge = 1;
            else
                baseFragment.iCharge = -1;
            vFragmentResults.push_back(baseFragment);
        }
    }
    return vFragmentResults;
}

bool OligoProcessUtils::ReadFastaAndDigest(ThreadPool *tp, vector<NucleicAcid>& vNucleicAcid, vector<Oligonucleotide>& vOligo)
{
   bool bSucceeded = true;

   bSucceeded = this->ReadFasta(Param::getInstance().options.sDatabaseName, vNucleicAcid);
   if (!bSucceeded) return false;

   bSucceeded = ModificationUtils::getInstance().AddStaticMod(vNucleicAcid, Param::getInstance().options.sStaticModsDatabase);
   if (!bSucceeded) return false;

   bSucceeded = this->ReadEnzymeFormula(Param::getInstance().options.sEnzymeDatabase, Param::getInstance().options.sEnzymeName);
   if (!bSucceeded) return false;

   for (int index = 0; index < vNucleicAcid.size(); index++)
   {
       NucleicAcid& nucleicAcid = vNucleicAcid[index];
      //  while (pSearchThreadPool->jobs_.size() >= CommonValues::iMaxThreadProc)
      //  {
      //      pSearchThreadPool->wait_for_available_thread();
      //  }
      //  pSearchThreadPool->doJob(std::bind(DigestNucleicAcid, nucleicAcid, std::ref(vOligonucleotide)));

      this->Digest(nucleicAcid, vOligo);
       
       if (!(index % 50))
       {
           //WriteOutUtils::ShowProgressBar(index, vNucleicAcid.size(), "Digest Nucleic Acid", "");
       }

   }
   // pSearchThreadPool->wait_on_threads();
   //WriteOutUtils::ShowProgressBar(vNucleicAcid.size(), vNucleicAcid.size(), "Digest Nucleic Acid", "done\n");

    if (!bSucceeded)
        return false;
    else
        return true;
}

bool OligoProcessUtils::ReadFasta(const std::string& fastaFilename, vector<NucleicAcid>& vNucleicAcid)
{
   ifstream fastaFile(fastaFilename);
   if (!fastaFile.is_open()) {
       cout << "Error - cannot read database file \"" + fastaFilename + "\"\n";
       return false;
   }

   string line;
   NucleicAcid nucleicAcid;
   int index = 0;
   while (getline(fastaFile, line)) {
       if (line.empty()) continue; // Skip empty lines
       if (line[0] == '>') { // Header line
           if (!nucleicAcid.sName.empty() && !nucleicAcid.sSequence.empty()) {
               vNucleicAcid.push_back(nucleicAcid);
               index++;
           }
           nucleicAcid.sName = line.substr(1); // Store the name without '>'
           nucleicAcid.sSequence.clear(); // Clear sequence for new entry
           nucleicAcid.iWhichNucleicAcid = index;
       } else {
           nucleicAcid.sSequence += line; // Append sequence lines
       }
   }
   vNucleicAcid.push_back(nucleicAcid); // Add the last nucleic acid read
   fastaFile.close();
   return true;
}

void OligoProcessUtils::Digest(const NucleicAcid& nucleicAcid, vector<Oligonucleotide>& vOligo) //会使用多线程
{
   if (this->vEnzymeFormula.empty())
   {
      std::cerr << "Error: Enzyme formula is not set" << std::endl;
      return;
   }

   //generate cut positions
   std::vector<int> vCutPositions;
   vCutPositions.push_back(0);
   for (size_t index = 1; index <= nucleicAcid.sSequence.size() - 1; index++)
   {
      for (const auto& formula : this->vEnzymeFormula)
      {
         if (formula.find('*') != std::string::npos)
         {
            std::string formulaLeftPart = formula.substr(0, formula.find('*')); // Get the left part of the formula
            std::string formulaRightPart = formula.substr(formula.find('*') + 1); // Get the right part of the formula
            std::string sequenceLeftPart = nucleicAcid.sSequence.substr(0, index); // Get the left part of the sequence
            std::string sequenceRightPart = nucleicAcid.sSequence.substr(index); // Get the right part of the sequence
            if (sequenceLeftPart.size() >= formulaLeftPart.size() && 
                sequenceLeftPart.substr(sequenceLeftPart.size() - formulaLeftPart.size()) == formulaLeftPart &&
                sequenceRightPart.size() >= formulaRightPart.size() && 
                sequenceRightPart.substr(0, formulaRightPart.size()) == formulaRightPart)
            {
               vCutPositions.push_back(index);
            }
         }
         else
            throw std::runtime_error("Enzyme formula is not valid");
      }
   }
   vCutPositions.push_back(nucleicAcid.sSequence.size());

   //generate oligonucleotides
   vector<Modification> vOligo3end = MassUtils::getInstance().GenerateEndModification(Param::getInstance().options.sOligo3Mod);
   vector<Modification> vOligo5end = MassUtils::getInstance().GenerateEndModification(Param::getInstance().options.sOligo5Mod);
   vector<Modification> vNucleicAcid3end = MassUtils::getInstance().GenerateEndModification(Param::getInstance().options.sNucleic3Mod); 
   vector<Modification> vNucleicAcid5end = MassUtils::getInstance().GenerateEndModification(Param::getInstance().options.sNucleic5Mod);

   for (size_t i = 0; i < vCutPositions.size() - 1; ++i) {
       int iStartPos = vCutPositions[i];
       for (size_t addIndex = 1; addIndex <= Param::getInstance().options.iAllowedMissed + 1; ++addIndex)
       {
         if (i + addIndex >= vCutPositions.size() - 1)
            break;
         int iEndPos = vCutPositions[i + addIndex];
         //check avaliable length
         if (iEndPos - iStartPos < Param::getInstance().options.oligoLengthRange.iStart || iEndPos - iStartPos > Param::getInstance().options.oligoLengthRange.iEnd)
            continue;
         
         std::string sSequence = nucleicAcid.sSequence.substr(iStartPos, iEndPos - iStartPos);
         if (this->mOligoIndex.find(sSequence) != this->mOligoIndex.end())
         {
            Oligonucleotide::NucleicAcidIndex whichNucleicAcid;
            whichNucleicAcid.iWhichNucleicAcid = nucleicAcid.iWhichNucleicAcid;
            whichNucleicAcid.iNABeginPos = iStartPos;
            whichNucleicAcid.iNAEndPos = iEndPos - 1;
            const std::vector<int>& vPositions = this->mOligoIndex[sSequence];
            for (int iPosition : vPositions)
            {
               if (iPosition >=0 && iPosition < vOligo.size() && vOligo[iPosition].sSequence == sSequence)
               {
                  vOligo[iPosition].vWhichNucleicAcid.push_back(whichNucleicAcid);
               }
            }
         }
         else{
            //generate candidate modification
            vector<Modification> vCandidate3end;
            vector<Modification> vCandidate5end;
            if (iStartPos == 0)
               vCandidate5end = vNucleicAcid5end;
            else
               vCandidate5end = vOligo5end;
            if (iEndPos == nucleicAcid.sSequence.size())
               vCandidate3end = vNucleicAcid3end;
            else
               vCandidate3end = vOligo3end;
            
            //generate oligonucleotide
            std::vector<int> vOligoIndex;
            for (auto& _5endMod : vCandidate5end)
            {
               for (auto& _3endMod : vCandidate3end)
               {
                  Oligonucleotide oligonucleotide;
                  Oligonucleotide::NucleicAcidIndex whichNucleicAcid;
                  whichNucleicAcid.iWhichNucleicAcid = nucleicAcid.iWhichNucleicAcid;
                  whichNucleicAcid.iNABeginPos = iStartPos;
                  whichNucleicAcid.iNAEndPos = iEndPos - 1;
                  oligonucleotide.sSequence = nucleicAcid.sSequence.substr(iStartPos, iEndPos - iStartPos);
                  oligonucleotide.vWhichNucleicAcid.push_back(whichNucleicAcid);
                  oligonucleotide.bDecoy = false;
                  oligonucleotide.end3TermMod = _3endMod;
                  oligonucleotide.end5TermMod = _5endMod;

                  for (const auto& mod : nucleicAcid.vMods) 
                  {
                     int modPos = mod.first;
                     if (modPos >= whichNucleicAcid.iNABeginPos && modPos <= whichNucleicAcid.iNAEndPos) {
                        oligonucleotide.vMods.push_back({modPos - whichNucleicAcid.iNABeginPos, mod.second});
                     }
                  }

                  OligonucleotideID oligoID;
                  double dCalcMass = MassUtils::getInstance().CalOligoMass(oligonucleotide, oligoID, true, true, false); //no consider of various modification (just for original mass check)
                  
                  if (Param::getInstance().options.digestMassRange.inRange(dCalcMass))
                  {
                     vOligoIndex.push_back(vOligo.size());
                     vOligo.push_back(oligonucleotide);
                  }
               }
            }
            this->mOligoIndex[sSequence] = vOligoIndex;
         }
       }
   }
   return;
}

bool OligoProcessUtils::ReadEnzymeFormula(const std::string& sEnzymeDatabase, const std::string& sEnzymeName)
{
   TSVParser tsvParser;
   if (!tsvParser.readFile(sEnzymeDatabase)) {
       std::cerr << "Error: Could not read enzyme database file " << sEnzymeDatabase << std::endl;
       return false;
   }

   std::vector<std::map<std::string, std::string>> data = tsvParser.getData();
   std::string sEnzymeFormula;
   for (const auto& row : data) {
      if (row.at("enzyme_name") == sEnzymeName) {
         sEnzymeFormula = row.at("enzyme_formula");
         break;
      }
   }
   
   size_t startPos = sEnzymeFormula.find('(');
   size_t endPos = sEnzymeFormula.find(')');
   if (startPos != std::string::npos && endPos != std::string::npos && endPos > startPos) {
       std::string formulaContent = sEnzymeFormula.substr(startPos + 1, endPos - startPos - 1);
       std::stringstream ss(formulaContent);
       std::string fragment;
       while (std::getline(ss, fragment, ',')) {
           fragment.erase(0, fragment.find_first_not_of(" \t")); // Remove leading whitespace
           fragment.erase(fragment.find_last_not_of(" \t") + 1); // Remove trailing whitespace
           this->vEnzymeFormula.push_back(fragment);
       }
   }
   else
   {
      std::cerr << "Error: Could not find enzyme formula for " << sEnzymeName << std::endl;
      return false;
   }

   return true;

}

bool OligoProcessUtils::GenerateDecoyOligo(vector<Oligonucleotide>& vOligo, int iDecoyOption)
{
   std::set<std::string> sDecoySequence;
   for (auto& oligo : vOligo)
      sDecoySequence.insert(oligo.sSequence);

   //generate decoy oligonucleotide (todo: shuffle may cause same decoy sequence)
    vector<Oligonucleotide> vPreserveList;
    for (auto& oligo : vOligo)
    {
        if (oligo.sSequence.size() <= 2)
            continue;
        Oligonucleotide decoy = oligo;
        decoy.bDecoy = true;
        decoy.vMods.clear();

         if (iDecoyOption == 0) //reverse
         {
            string sReverseSequence = oligo.sSequence.substr(0, oligo.sSequence.size() - 1);
            reverse(sReverseSequence.begin(), sReverseSequence.end());
            decoy.sSequence = sReverseSequence + oligo.sSequence.back();

            for (auto& mod : oligo.vMods)
            {
               int index;
               // 如果原本是最后一个字符的 modification，不进行位置反转
               if (mod.first == oligo.sSequence.size() - 1)
                  index = mod.first;  // 保持最后一个字符的位置不变
               else
                  // 对于其他字符，调整位置
                  index = oligo.sSequence.size() - 1 - mod.first;
               decoy.vMods.push_back({ index, mod.second });
            }

            if (sDecoySequence.find(decoy.sSequence) == sDecoySequence.end())
            vPreserveList.push_back(decoy);
         }
         else if (iDecoyOption == 1) //shuffle
         {
            for (size_t repeat_index = 0; repeat_index < 4; repeat_index++)
            {
               int iSeqSize = oligo.sSequence.size();
               vector<int> vShuffleIndices(iSeqSize - 1);
               for (int index = 0; index < vShuffleIndices.size(); index++) // 除去最后一个字符
               {
                  vShuffleIndices[index] = index;
               }

               // 随机打乱序列的索引
               random_device rd;
               mt19937 g(rd());
               shuffle(vShuffleIndices.begin(), vShuffleIndices.end(), g);

               string reverseSequence;
               for (int& index : vShuffleIndices)
               {
                  reverseSequence.push_back(oligo.sSequence[index]);
               }
               reverseSequence.push_back(oligo.sSequence.back());
               decoy.sSequence = reverseSequence;

               // 更新 modifications 的位置
               for (const auto& mod : oligo.vMods)
               {
                  int originalPos = mod.first;
                  // 如果修改位置是最后一个字符，直接保持不变
                  if (originalPos == oligo.sSequence.size() - 1) {
                     decoy.vMods.push_back({ originalPos, mod.second });
                  }
                  else {
                     // 如果是其他字符，根据 shuffle 后的索引调整位置
                     auto it = find(vShuffleIndices.begin(), vShuffleIndices.end(), originalPos);
                     if (it != vShuffleIndices.end()) {
                           int newPos = distance(vShuffleIndices.begin(), it);
                           decoy.vMods.push_back({ newPos, mod.second });
                     }
                  }
               }

               if (sDecoySequence.find(decoy.sSequence) == sDecoySequence.end())
               {
                  vPreserveList.push_back(decoy);
               }
            }
         }
    }

    for (Oligonucleotide& oligo : vPreserveList)
    {
        vOligo.push_back(oligo);
    }

    return true;

}
