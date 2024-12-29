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

void MassUtils::AssignMass(MassUtil& massUtility)
{
   double H, O, C, N, P, S, Se;
   if (massUtility.bMonoMass) // monoisotopic masses
   {
      H = massUtility.pdAtomMass[(int)'h'] = ElementMassUtils::Hydrogen_Mono; // hydrogen
      O = massUtility.pdAtomMass[(int)'o'] = ElementMassUtils::Oxygen_Mono;  // oxygen
      C = massUtility.pdAtomMass[(int)'c'] = ElementMassUtils::Carbon_Mono;   // carbon
      N = massUtility.pdAtomMass[(int)'n'] = ElementMassUtils::Nitrogen_Mono;   // nitrogen
      P = massUtility.pdAtomMass[(int)'p'] = ElementMassUtils::Phosphorus_Mono;    // phosphorus
      S = massUtility.pdAtomMass[(int)'s'] = ElementMassUtils::Sulfur_Mono;   // sulphur
      Se = massUtility.pdAtomMass[(int)'e'] = ElementMassUtils::Selenium_Mono;  // selenium
   }
   else  // average masses
   {
      H = massUtility.pdAtomMass[(int)'h'] = ElementMassUtils::Hydrogen_Ave;
      O = massUtility.pdAtomMass[(int)'o'] = ElementMassUtils::Oxygen_Ave;
      C = massUtility.pdAtomMass[(int)'c'] = ElementMassUtils::Carbon_Ave;
      N = massUtility.pdAtomMass[(int)'n'] = ElementMassUtils::Nitrogen_Ave;
      P = massUtility.pdAtomMass[(int)'p'] = ElementMassUtils::Phosphorus_Ave;
      S = massUtility.pdAtomMass[(int)'s'] = ElementMassUtils::Sulfur_Ave;
      Se = massUtility.pdAtomMass[(int)'e'] = ElementMassUtils::Selenium_Ave;
   }

   double backbone;
   if (massUtility.bDeoxidation)
   {
      backbone = massUtility.dBackbone = C*5 + H*8 + O*5 + P;
   }
   else
   {
      backbone = massUtility.dBackbone = C*5 + H*8 + O*6 + P;
   }

   massUtility.pdBaseMass[(int)'A'] = C*5  + H*4  + N*5;
   massUtility.pdBaseMass[(int)'G'] = C*5  + H*4  + N*5  + O;
   massUtility.pdBaseMass[(int)'C'] = C*4  + H*4  + N*3  + O;
   massUtility.pdBaseMass[(int)'U'] = C*4  + H*3  + N*2  + O*2;
   massUtility.pdBaseMass[(int)'T'] = C*5  + H*5  + N*2  + O*2;

   massUtility.pdNTMass[(int)'A'] = massUtility.pdBaseMass[(int)'A'] + backbone;
   massUtility.pdNTMass[(int)'G'] = massUtility.pdBaseMass[(int)'G'] + backbone;
   massUtility.pdNTMass[(int)'C'] = massUtility.pdBaseMass[(int)'C'] + backbone;
   massUtility.pdNTMass[(int)'U'] = massUtility.pdBaseMass[(int)'U'] + backbone;
   massUtility.pdNTMass[(int)'T'] = massUtility.pdBaseMass[(int)'T'] + backbone;
}

void MassUtils::AssignAAMass(MassUtil& massUtility)
{
   double H, O, C, N, S, Se;

   if (massUtility.bMonoMass) // monoisotopic masses
   {
      H = ElementMassUtils::Hydrogen_Mono; // hydrogen
      O = ElementMassUtils::Oxygen_Mono;  // oxygen
      C = ElementMassUtils::Carbon_Mono;   // carbon
      N = ElementMassUtils::Nitrogen_Mono;   // nitrogen
      S = ElementMassUtils::Sulfur_Mono;   // sulphur
      Se = ElementMassUtils::Selenium_Mono;  // selenium
   }
   else  // average masses
   {
      H = ElementMassUtils::Hydrogen_Ave; // hydrogen
      O = ElementMassUtils::Oxygen_Ave;  // oxygen
      C = ElementMassUtils::Carbon_Ave;   // carbon
      N = ElementMassUtils::Nitrogen_Ave;   // nitrogen
      S = ElementMassUtils::Sulfur_Ave;   // sulphur
      Se = ElementMassUtils::Selenium_Ave;  // selenium
   }

   massUtility.pdNTMass[(int)'G'] = C*2  + H*3  + N   + O ;
   massUtility.pdNTMass[(int)'A'] = C*3  + H*5  + N   + O ;
   massUtility.pdNTMass[(int)'S'] = C*3  + H*5  + N   + O*2 ;
   massUtility.pdNTMass[(int)'P'] = C*5  + H*7  + N   + O ;
   massUtility.pdNTMass[(int)'V'] = C*5  + H*9  + N   + O ;
   massUtility.pdNTMass[(int)'T'] = C*4  + H*7  + N   + O*2 ;
   massUtility.pdNTMass[(int)'C'] = C*3  + H*5  + N   + O   + S ;
   massUtility.pdNTMass[(int)'U'] = C*3  + H*5  + N   + O   + Se ;
   massUtility.pdNTMass[(int)'L'] = C*6  + H*11 + N   + O ;
   massUtility.pdNTMass[(int)'I'] = C*6  + H*11 + N   + O ;
   massUtility.pdNTMass[(int)'N'] = C*4  + H*6  + N*2 + O*2 ;
   massUtility.pdNTMass[(int)'D'] = C*4  + H*5  + N   + O*3 ;
   massUtility.pdNTMass[(int)'Q'] = C*5  + H*8  + N*2 + O*2 ;
   massUtility.pdNTMass[(int)'K'] = C*6  + H*12 + N*2 + O ;
   massUtility.pdNTMass[(int)'E'] = C*5  + H*7  + N   + O*3 ;
   massUtility.pdNTMass[(int)'M'] = C*5  + H*9  + N   + O   + S ;
   massUtility.pdNTMass[(int)'H'] = C*6  + H*7  + N*3 + O ;
   massUtility.pdNTMass[(int)'F'] = C*9  + H*9  + N   + O ;
   massUtility.pdNTMass[(int)'R'] = C*6  + H*12 + N*4 + O ;
   massUtility.pdNTMass[(int)'Y'] = C*9  + H*9  + N   + O*2 ;
   massUtility.pdNTMass[(int)'W'] = C*11 + H*10 + N*2 + O ;
   massUtility.pdNTMass[(int)'O'] = C*12  + H*19 + N*3 + O*2 ;
}

double MassUtils::CalculateOligonucleotideMass(Oligonucleotide& oligonucleotide, bool bAddTermMod, bool bAddStaticMod, int iWhichVarModCombination)
{
   double dCalcMass = 0.0;

   if (oligonucleotide.sSequence.size() == 0)
      return dCalcMass;


   for (char& c : oligonucleotide.sSequence)
   {
      dCalcMass += (double) Param::g_staticParams.massUtility.pdNTMass[(int)c];
   }

   if (bAddTermMod)
   {
      dCalcMass += oligonucleotide.end3TermMod.dModMass;
      dCalcMass += oligonucleotide.end5TermMod.dModMass;
   }

   if (bAddStaticMod)
   {
      for (auto& staticModification : oligonucleotide.vMods)
      {
         dCalcMass += staticModification.second.dModMass;
      }
   }

   if (iWhichVarModCombination >= 0)
   {
      dCalcMass += Param::g_vVarModification[oligonucleotide.iVarModificationIndex].vVarModificationMassList[iWhichVarModCombination];
   }

   oligonucleotide.dMass = dCalcMass;
   return dCalcMass;

}

double MassUtils::CalculateOligonucleotideFragmentMass(OligonucleotideFragment& fragment, bool bAddTermMod, bool bAddMod)
{
   double dCalcMass = 0.0;

   if (fragment.sSequence.size() == 0)
      return dCalcMass;


   for (char& c : fragment.sSequence)
   {
      dCalcMass += (double) Param::g_staticParams.massUtility.pdNTMass[(int)c];
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
   
   fragment.dMass = dCalcMass;
   return dCalcMass;
}

double MassUtils::MassToMZ(double dMass, int iCharge)
{
   double dMZ = 0.0;
   double dAdduct = 0.0;

   if (Param::g_staticParams.options.sAdductIon == "H-" && !Param::g_staticParams.options.bPositive)
       dAdduct = ElementMassUtils::Proton_Mono;
   else if (Param::g_staticParams.options.sAdductIon == "H+" && Param::g_staticParams.options.bPositive)
       dAdduct = ElementMassUtils::Proton_Mono;
   else if (Param::g_staticParams.options.sAdductIon == "K+" && Param::g_staticParams.options.bPositive)
       dAdduct = ElementMassUtils::Potassium_Mono - (ElementMassUtils::Hydrogen_Mono - ElementMassUtils::Proton_Mono);
   else if (Param::g_staticParams.options.sAdductIon == "Na+" && Param::g_staticParams.options.bPositive)
       dAdduct = ElementMassUtils::Sodium_Mono - (ElementMassUtils::Hydrogen_Mono - ElementMassUtils::Proton_Mono);
   else if (Param::g_staticParams.options.sAdductIon == "NH4+" && Param::g_staticParams.options.bPositive)
       dAdduct = (ElementMassUtils::Nitrogen_Mono + 3 * ElementMassUtils::Hydrogen_Mono + ElementMassUtils::Proton_Mono);
   else
   {
       string sErrorMsg = "Error: not valid adduct ion selected, may cause error in calculating mz;";
       logerr(sErrorMsg.c_str());
       dAdduct = 0.0;
   }

   if (iCharge == 0)
      dMZ = dMass;
   else
      dMZ = (dMass + iCharge * dAdduct) / fabs((double)iCharge);
   
   if (dMZ < 0)
      dMZ = 0.0;

   return dMZ;
}

double MassUtils::MZToMass(double dMZ, int iCharge)
{
   double dMass = 0.0;
   double dAdduct = 0.0;

   if (Param::g_staticParams.options.sAdductIon == "H-" && !Param::g_staticParams.options.bPositive)
       dAdduct = ElementMassUtils::Proton_Mono;
   else if (Param::g_staticParams.options.sAdductIon == "H+" && Param::g_staticParams.options.bPositive)
       dAdduct = ElementMassUtils::Proton_Mono;
   else if (Param::g_staticParams.options.sAdductIon == "K+" && Param::g_staticParams.options.bPositive)
       dAdduct = ElementMassUtils::Potassium_Mono - (ElementMassUtils::Hydrogen_Mono - ElementMassUtils::Proton_Mono);
   else if (Param::g_staticParams.options.sAdductIon == "Na+" && Param::g_staticParams.options.bPositive)
       dAdduct = ElementMassUtils::Sodium_Mono - (ElementMassUtils::Hydrogen_Mono - ElementMassUtils::Proton_Mono);
   else if (Param::g_staticParams.options.sAdductIon == "NH4+" && Param::g_staticParams.options.bPositive)
       dAdduct = (ElementMassUtils::Nitrogen_Mono + 3 * ElementMassUtils::Hydrogen_Mono + ElementMassUtils::Proton_Mono);
   else
   {
       string sErrorMsg = "Error: not valid adduct ion selected, may cause error in calculating mz;";
       logerr(sErrorMsg.c_str());
       dAdduct = 0.0;
   }

   if (iCharge == 0)
      dMass = dMZ;
   else
      dMass = dMZ * (fabs((double)iCharge)) - iCharge * dAdduct;
   
   if (dMass < 0)
      dMass = 0.0;

   return dMass;
}
