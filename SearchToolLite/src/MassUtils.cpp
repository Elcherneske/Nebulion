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

#include "MassUtils.h"
#include "Param.h"

#include <string>
#include <sstream>


MassUtils::MassUtils()
{
   double H, O, C, N, P, S, Se;
   H = ElementMassUtils::Hydrogen_Mono; // hydrogen
   O = ElementMassUtils::Oxygen_Mono;  // oxygen
   C = ElementMassUtils::Carbon_Mono;   // carbon
   N = ElementMassUtils::Nitrogen_Mono;   // nitrogen
   P = ElementMassUtils::Phosphorus_Mono;    // phosphorus
   S = ElementMassUtils::Sulfur_Mono;   // sulphur
   Se = ElementMassUtils::Selenium_Mono;  // selenium

   // nucleotide backbone mass
   if (Param::getInstance().options.bDeoxidation)
      dBackbone = C*5 + H*8 + O*5 + P;
   else
      dBackbone = C*5 + H*8 + O*6 + P;

   // nucleotide base mass
   std::fill(std::begin(pdBaseMass), std::end(pdBaseMass), 99999);
   std::fill(std::begin(pdNTMass), std::end(pdNTMass), 99999);
   
   pdBaseMass[(int)'A'] = C*5  + H*4  + N*5;
   pdBaseMass[(int)'G'] = C*5  + H*4  + N*5  + O;
   pdBaseMass[(int)'C'] = C*4  + H*4  + N*3  + O;
   pdBaseMass[(int)'U'] = C*4  + H*3  + N*2  + O*2;
   pdBaseMass[(int)'T'] = C*5  + H*5  + N*2  + O*2;

   pdNTMass[(int)'A'] = pdBaseMass[(int)'A'] + dBackbone;
   pdNTMass[(int)'G'] = pdBaseMass[(int)'G'] + dBackbone;
   pdNTMass[(int)'C'] = pdBaseMass[(int)'C'] + dBackbone;
   pdNTMass[(int)'U'] = pdBaseMass[(int)'U'] + dBackbone;
   pdNTMass[(int)'T'] = pdBaseMass[(int)'T'] + dBackbone;

   // adduct mass
   mapAdductMass["H+"] = ElementMassUtils::Proton_Mono;
   mapAdductMass["K+"] = ElementMassUtils::Potassium_Mono - (ElementMassUtils::Hydrogen_Mono - ElementMassUtils::Proton_Mono);
   mapAdductMass["Na+"] = ElementMassUtils::Sodium_Mono - (ElementMassUtils::Hydrogen_Mono - ElementMassUtils::Proton_Mono);
   mapAdductMass["NH4+"] = (ElementMassUtils::Nitrogen_Mono + 3 * ElementMassUtils::Hydrogen_Mono + ElementMassUtils::Proton_Mono);
}

MassUtils::~MassUtils()
{
}


// void MassUtils::AssignAAMass(MassUtil& massUtility)
// {
//    double H, O, C, N, S, Se;

//    if (massUtility.bMonoMass) // monoisotopic masses
//    {
//       H = ElementMassUtils::Hydrogen_Mono; // hydrogen
//       O = ElementMassUtils::Oxygen_Mono;  // oxygen
//       C = ElementMassUtils::Carbon_Mono;   // carbon
//       N = ElementMassUtils::Nitrogen_Mono;   // nitrogen
//       S = ElementMassUtils::Sulfur_Mono;   // sulphur
//       Se = ElementMassUtils::Selenium_Mono;  // selenium
//    }
//    else  // average masses
//    {
//       H = ElementMassUtils::Hydrogen_Ave; // hydrogen
//       O = ElementMassUtils::Oxygen_Ave;  // oxygen
//       C = ElementMassUtils::Carbon_Ave;   // carbon
//       N = ElementMassUtils::Nitrogen_Ave;   // nitrogen
//       S = ElementMassUtils::Sulfur_Ave;   // sulphur
//       Se = ElementMassUtils::Selenium_Ave;  // selenium
//    }

//    massUtility.pdNTMass[(int)'G'] = C*2  + H*3  + N   + O ;
//    massUtility.pdNTMass[(int)'A'] = C*3  + H*5  + N   + O ;
//    massUtility.pdNTMass[(int)'S'] = C*3  + H*5  + N   + O*2 ;
//    massUtility.pdNTMass[(int)'P'] = C*5  + H*7  + N   + O ;
//    massUtility.pdNTMass[(int)'V'] = C*5  + H*9  + N   + O ;
//    massUtility.pdNTMass[(int)'T'] = C*4  + H*7  + N   + O*2 ;
//    massUtility.pdNTMass[(int)'C'] = C*3  + H*5  + N   + O   + S ;
//    massUtility.pdNTMass[(int)'U'] = C*3  + H*5  + N   + O   + Se ;
//    massUtility.pdNTMass[(int)'L'] = C*6  + H*11 + N   + O ;
//    massUtility.pdNTMass[(int)'I'] = C*6  + H*11 + N   + O ;
//    massUtility.pdNTMass[(int)'N'] = C*4  + H*6  + N*2 + O*2 ;
//    massUtility.pdNTMass[(int)'D'] = C*4  + H*5  + N   + O*3 ;
//    massUtility.pdNTMass[(int)'Q'] = C*5  + H*8  + N*2 + O*2 ;
//    massUtility.pdNTMass[(int)'K'] = C*6  + H*12 + N*2 + O ;
//    massUtility.pdNTMass[(int)'E'] = C*5  + H*7  + N   + O*3 ;
//    massUtility.pdNTMass[(int)'M'] = C*5  + H*9  + N   + O   + S ;
//    massUtility.pdNTMass[(int)'H'] = C*6  + H*7  + N*3 + O ;
//    massUtility.pdNTMass[(int)'F'] = C*9  + H*9  + N   + O ;
//    massUtility.pdNTMass[(int)'R'] = C*6  + H*12 + N*4 + O ;
//    massUtility.pdNTMass[(int)'Y'] = C*9  + H*9  + N   + O*2 ;
//    massUtility.pdNTMass[(int)'W'] = C*11 + H*10 + N*2 + O ;
//    massUtility.pdNTMass[(int)'O'] = C*12  + H*19 + N*3 + O*2 ;
// }

double MassUtils::CalOligoMass(const Oligonucleotide& oligo, const OligonucleotideID& oligoID, bool bAddTermMod, bool bAddStaticMod, bool bAddVarMod)
{
   double dCalcMass = 0.0;
   if (oligo.sSequence.size() == 0)
      return dCalcMass;

   for (const char& c : oligo.sSequence)
   {
      dCalcMass += pdNTMass[(int)c];
   }

   if (bAddTermMod)
   {
      dCalcMass += oligo.end3TermMod.dModMass;
      dCalcMass += oligo.end5TermMod.dModMass;
   }

   if (bAddStaticMod)
   {
      for (const auto& mod : oligo.vMods)
      {
         dCalcMass += mod.second.dModMass;
      }
   }

   if (bAddVarMod)
   {
      dCalcMass += Param::getInstance().vVarModFormulas[oligoID.iVarModFormulaIndex].vVarModMassList[oligoID.iVarModIndex];
   }

   return dCalcMass;

}

double MassUtils::CalFragmentMass(const OligoFragment& fragment, bool bAddTermMod, bool bAddMod)
{
   double dCalcMass = 0.0;

   if (fragment.sSequence.size() == 0)
      return dCalcMass;


   for (const char& c : fragment.sSequence)
   {
      dCalcMass += pdNTMass[(int)c];
   }

   if (bAddTermMod)
   {
      dCalcMass += fragment.end3TermMod.dModMass;
      dCalcMass += fragment.end5TermMod.dModMass;
   }

   if (bAddMod)
   {
      for (auto& staticModification : fragment.vMods)
      {
         dCalcMass += staticModification.second.dModMass;
      }
   }
   return dCalcMass;
}

std::vector<Modification> MassUtils::GenerateEndModification(std::string sCandidates)
{
   std::vector<Modification> vModification;
   std::vector<std::string> vCandidates;
   std::stringstream ss(sCandidates);
   std::string sCandidate;
   while (std::getline(ss, sCandidate, ','))
   {
      sCandidate.erase(0, sCandidate.find_first_not_of(" \t")); // Remove leading whitespace
      sCandidate.erase(sCandidate.find_last_not_of(" \t") + 1); // Remove trailing whitespace
      vCandidates.push_back(sCandidate);
   }

   for (std::string sCandidate : vCandidates)
   {
      Modification modification;
      modification.sModName = sCandidate;
      bool bDelete;
      if (sCandidate.front() == '@')
      {
         bDelete = true;
         sCandidate = sCandidate.substr(1);
      }
      else
         bDelete = false;
      double dModMass = 0.0;
      if (sCandidate == "H")
         dModMass = ElementMassUtils::Hydrogen_Mono;
      else if (sCandidate == "OH")
         dModMass = ElementMassUtils::Hydrogen_Mono + ElementMassUtils::Oxygen_Mono;
      else if (sCandidate == "H2PO3")
         dModMass = ElementMassUtils::Hydrogen_Mono * 2 + ElementMassUtils::Phosphorus_Mono + ElementMassUtils::Oxygen_Mono * 3;
      else if (sCandidate == "NH2")
         dModMass = ElementMassUtils::Nitrogen_Mono + ElementMassUtils::Hydrogen_Mono * 2;
      else if (sCandidate == "PO2")
          dModMass = ElementMassUtils::Phosphorus_Mono + ElementMassUtils::Oxygen_Mono * 2;
      else
         dModMass = 0.0;
      
      if (bDelete)
         dModMass = -dModMass;
         
      modification.dModMass = dModMass;
      vModification.push_back(modification);
   }

   return vModification;
}

double MassUtils::CalCutModMass(IonSpecies& ionSpecies)
{
   double dModMass = 0.0;
   double C = ElementMassUtils::Carbon_Mono;
   double H = ElementMassUtils::Hydrogen_Mono;
   double O = ElementMassUtils::Oxygen_Mono;
   double P = ElementMassUtils::Phosphorus_Mono;
   switch (ionSpecies.ionType)
   {
      case IonSpecies::ion_a:
         dModMass -= (H + 3 * O + P);
         if (!ionSpecies.bFreeRadical)
            dModMass -= H;
         return dModMass;

      case IonSpecies::ion_b:
         dModMass -= (H + 2 * O + P);
         if (!ionSpecies.bFreeRadical)
            dModMass += H;
         return dModMass;

      case IonSpecies::ion_c:
         dModMass -= 0.0;
         if (!ionSpecies.bFreeRadical)
            dModMass -= H;
         return dModMass;

      case IonSpecies::ion_d:
         dModMass += O;
         if (!ionSpecies.bFreeRadical)
            dModMass += H;
         return dModMass;
      
      case IonSpecies::ion_w:
         dModMass += (H + 3 * O + P);
         if (!ionSpecies.bFreeRadical)
            dModMass += H;
         return dModMass;
      
      case IonSpecies::ion_x:
         dModMass += (H + 2 * O + P);
         if (!ionSpecies.bFreeRadical)
            dModMass -= H;
         return dModMass;

      case IonSpecies::ion_y:
         dModMass -= 0.0;
         if (!ionSpecies.bFreeRadical)
            dModMass += H;
         return dModMass;

      case IonSpecies::ion_z:
         dModMass -= O;
         if (!ionSpecies.bFreeRadical)
            dModMass -= H;
         return dModMass;
      default: //NL_ion
         return 0;
   }
}

double MassUtils::CalCutModMass(IonSpecies& ionSpecies, char Base)
{
   double dModMass = 0.0;
   double C = ElementMassUtils::Carbon_Mono;
   double H = ElementMassUtils::Hydrogen_Mono;
   double O = ElementMassUtils::Oxygen_Mono;
   double P = ElementMassUtils::Phosphorus_Mono;
   switch (ionSpecies.ionType)
   {
      case IonSpecies::ion_a:
         dModMass -= (H + 3 * O + P);
         if (ionSpecies.bDropBase)
            dModMass -= pdBaseMass[(int)Base];
         if (!ionSpecies.bFreeRadical)
            dModMass -= H;
         return dModMass;

      case IonSpecies::ion_b:
         dModMass -= (H + 2 * O + P);
         if (ionSpecies.bDropBase)
            dModMass -= pdBaseMass[(int)Base];
         if (!ionSpecies.bFreeRadical)
            dModMass += H;
         return dModMass;

      case IonSpecies::ion_c:
         dModMass -= 0.0;
         if (ionSpecies.bDropBase)
            dModMass -= pdBaseMass[(int)Base];
         if (!ionSpecies.bFreeRadical)
            dModMass -= H;
         return dModMass;

      case IonSpecies::ion_d:
         dModMass += O;
         if (ionSpecies.bDropBase)
            dModMass -= pdBaseMass[(int)Base];
         if (!ionSpecies.bFreeRadical)
            dModMass += H;
         return dModMass;
      
      case IonSpecies::ion_w:
         dModMass += (H + 3 * O + P);
         if (ionSpecies.bDropBase)
            dModMass -= pdBaseMass[(int)Base];
         if (!ionSpecies.bFreeRadical)
            dModMass += H;
         return dModMass;
      
      case IonSpecies::ion_x:
         dModMass += (H + 2 * O + P);
         if (ionSpecies.bDropBase)
            dModMass -= pdBaseMass[(int)Base];
         if (!ionSpecies.bFreeRadical)
            dModMass -= H;
         return dModMass;

      case IonSpecies::ion_y:
         dModMass -= 0.0;
         if (ionSpecies.bDropBase)
            dModMass -= pdBaseMass[(int)Base];
         if (!ionSpecies.bFreeRadical)
            dModMass += H;
         return dModMass;

      case IonSpecies::ion_z:
         dModMass -= O;
         if (ionSpecies.bDropBase)
            dModMass -= pdBaseMass[(int)Base];
         if (!ionSpecies.bFreeRadical)
            dModMass -= H;
         return dModMass;
      default: //NL_ion
         return 0;
   }
}

double MassUtils::MassToMZ(double dMass, int iCharge, std::string sAdductIon)
{
   bool bDelete;
   if (sAdductIon.front() == '@') {
      bDelete = true;
      sAdductIon = sAdductIon.substr(1);
   } else
      bDelete = false;
   
   double dAdductMass = 0.0;
   auto it = mapAdductMass.find(sAdductIon);
   if (it != mapAdductMass.end())
      dAdductMass = it->second;
   else
      throw std::invalid_argument("Adduct ion not found in map.");
   
   double dMZ = 0.0;
   if (iCharge == 0)
      dMZ = dMass;
   else
      dMZ = bDelete? (dMass - fabs((double)iCharge) * dAdductMass) / fabs((double)iCharge): (dMass + fabs((double)iCharge) * dAdductMass) / fabs((double)iCharge);
   
   if (dMZ < 0)
      throw std::invalid_argument("Calculated MZ cannot be negative.");

   return dMZ;
}

double MassUtils::MZToMass(double dMZ, int iCharge, std::string sAdductIon)
{
   bool bDelete;
   if (sAdductIon.front() == '@') {
      bDelete = true;
      sAdductIon = sAdductIon.substr(1);
   } else
      bDelete = false;
   
   double dAdductMass = 0.0;
   auto it = mapAdductMass.find(sAdductIon);
   if (it != mapAdductMass.end())
      dAdductMass = it->second;
   else
      throw std::invalid_argument("Adduct ion not found in map.");
   
   double dMass = 0.0;
   if (iCharge == 0)
      dMass = dMZ;
   else
      dMass = bDelete? (dMZ * (fabs((double)iCharge)) + fabs((double)iCharge) * dAdductMass) : (dMZ * (fabs((double)iCharge)) - fabs((double)iCharge) * dAdductMass);
   
   if (dMass < 0)
      throw std::invalid_argument("Calculated mass cannot be negative.");

   return dMass;
}

double MassUtils::GetBaseMass(char cBase)
{
   return pdBaseMass[(int)cBase];
}

double MassUtils::GetNTMass(char cBase)
{
   return pdNTMass[(int)cBase];
}

int MassUtils::CalBinIndex(double dMZ)
{
   return static_cast<int>(dMZ / Param::getInstance().options.dBinWidth);
}
