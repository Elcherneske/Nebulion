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

unordered_map<string, int>              OligoProcessUtils::decoyDuplicateMap;

OligoProcessUtils::OligoProcessUtils()
{
}

OligoProcessUtils::~OligoProcessUtils()
{
}

vector<OligonucleotideFragment> OligoProcessUtils::GenerateOligonucleotideFragment(const Oligonucleotide& oligonucleotide, IonSeries::Ions ionType, int iCharge, vector<pair<int, Modification>> vVarModList)
{
    vector<OligonucleotideFragment> vFragmentResults;
   if (ionType == IonSeries::ion_NL)
      return OligoProcessUtils::GenerateNLOligonucleotideFragment(oligonucleotide, iCharge);

   int iLenFragment = 0;
   int iLenOligonucleotide = oligonucleotide.sSequence.size();
   int iPosReverse;                       // points to residue in reverse order

   Modification bIonMod;
   Modification yIonMod;
   if (IonSeries::getBYTerm(ionType))
      bIonMod.dModMass = IonSeries::getTermMassMod(ionType);
   else
      yIonMod.dModMass = IonSeries::getTermMassMod(ionType);

   for (int i = 2; i < oligonucleotide.sSequence.size(); ++i)
   {

      iPosReverse = oligonucleotide.sSequence.size() - i;

      int IonLen = i;

      for (int charge = 1; charge <= abs(iCharge); charge++)
      {
         OligonucleotideFragment bIon;
         OligonucleotideFragment yIon;

         bIon.sSequence = oligonucleotide.sSequence.substr(0, IonLen);
         bIon.end5TermMod = oligonucleotide.end5TermMod;
         bIon.end3TermMod = bIonMod;
         if (ionType == IonSeries::ion_a_B)
         {
            bIon.end3TermMod.dModMass -= Param::g_staticParams.massUtility.pdBaseMass[(int)bIon.sSequence[bIon.sSequence.size() - 1]];
            bIon.end3TermMod.dModMass -= ElementMassUtils::Hydrogen_Mono;
         }
         bIon.ionType = ionType;
         if (iCharge > 0)
            bIon.iCharge = charge;
         else
            bIon.iCharge = -1 * charge;

         yIon.sSequence = oligonucleotide.sSequence.substr(iPosReverse, IonLen);
         yIon.end5TermMod = yIonMod;
         yIon.end3TermMod = oligonucleotide.end3TermMod;
         yIon.ionType = ionType;
         if (Param::g_staticParams.options.bPositive)
            yIon.iCharge = charge;
         else
            yIon.iCharge = -1 * charge;

         for (const auto& mod : oligonucleotide.vMods) {
            int modPos = mod.first;
            if (modPos >= 0 && modPos < i)
               bIon.vMods.push_back({modPos, mod.second});

            if (modPos >= iPosReverse && modPos < oligonucleotide.sSequence.size())
               yIon.vMods.push_back({modPos - iPosReverse, mod.second});
         }

         for (const auto& mod : vVarModList) {
            int modPos = mod.first;
            if (modPos >= 0 && modPos < i)
               bIon.vMods.push_back({modPos, mod.second});

            if (modPos >= iPosReverse && modPos < oligonucleotide.sSequence.size())
               yIon.vMods.push_back({modPos - iPosReverse, mod.second});
         }

         //todo: need a thread lock
         if (IonSeries::getBYTerm(ionType))
         {
            MassUtils::CalculateOligonucleotideFragmentMass(bIon, true, true);
            vFragmentResults.push_back(bIon);
         }
         else
         {
            MassUtils::CalculateOligonucleotideFragmentMass(yIon, true, true);
            vFragmentResults.push_back(yIon);
         }

      }
   }
   return vFragmentResults;
}

vector<OligonucleotideFragment> OligoProcessUtils::GenerateNLOligonucleotideFragment(const Oligonucleotide& oligonucleotide, int iCharge)
{
    vector<OligonucleotideFragment> vFragmentResults;

    double H = ElementMassUtils::Hydrogen_Mono;
    double P = ElementMassUtils::Phosphorus_Mono;
    double O = ElementMassUtils::Oxygen_Mono;
    OligonucleotideFragment fragment;
    fragment.end3TermMod = oligonucleotide.end3TermMod;
    fragment.end5TermMod = oligonucleotide.end5TermMod;
    fragment.ionType = IonSeries::ion_NL;
    fragment.vMods = oligonucleotide.vMods;

    // loss P or H2O
    for(int i = 1; i <= abs(iCharge); i++)
    {
        fragment.sSequence = oligonucleotide.sSequence;
        fragment.iCharge = i * iCharge / abs(iCharge);

        fragment.dMass = oligonucleotide.dMass - (H + P + 3 * O);
        vFragmentResults.push_back(fragment);

        fragment.dMass = oligonucleotide.dMass - (2 * H + O);
        vFragmentResults.push_back(fragment);

        fragment.dMass = oligonucleotide.dMass - (3 * H + P + 4 * O);
        vFragmentResults.push_back(fragment);
    }

    // base fragment and M-Base fragment
    if (oligonucleotide.sSequence.find('A') != string::npos)
    {
        //base fragment
        fragment.sSequence = "A";
        fragment.dMass = Param::g_staticParams.massUtility.pdBaseMass[(int)'A'] + ElementMassUtils::Hydrogen_Mono;
        if (Param::g_staticParams.options.bPositive)
            fragment.iCharge = 1;
        else 
            fragment.iCharge = -1;
        vFragmentResults.push_back(fragment);

        //M-base fragment
        for (int i = 1; i <= abs(iCharge); i++)
        {
            fragment.sSequence = oligonucleotide.sSequence;
            fragment.iCharge = i * iCharge / abs(iCharge);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'A'] + ElementMassUtils::Hydrogen_Mono);
            vFragmentResults.push_back(fragment);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'A'] + ElementMassUtils::Hydrogen_Mono) - (H + P + 3 * O);
            vFragmentResults.push_back(fragment);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'A'] + ElementMassUtils::Hydrogen_Mono) - (2 * H + O);
            vFragmentResults.push_back(fragment);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'A'] + ElementMassUtils::Hydrogen_Mono) - (3 * H + P + 4 * O);
            vFragmentResults.push_back(fragment);
        }
    }
    if (oligonucleotide.sSequence.find('T') != string::npos)
    {
        fragment.sSequence = "T";
        fragment.dMass = Param::g_staticParams.massUtility.pdBaseMass[(int)'T'] + ElementMassUtils::Hydrogen_Mono;
        if (Param::g_staticParams.options.bPositive)
            fragment.iCharge = 1;
        else 
            fragment.iCharge = -1;
      
        //M-base fragment
        for (int i = 1; i <= abs(iCharge); i++)
        {
            fragment.sSequence = oligonucleotide.sSequence;
            fragment.iCharge = i * iCharge / abs(iCharge);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'T'] + ElementMassUtils::Hydrogen_Mono);
            vFragmentResults.push_back(fragment);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'T'] + ElementMassUtils::Hydrogen_Mono) - (H + P + 3 * O);
            vFragmentResults.push_back(fragment);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'T'] + ElementMassUtils::Hydrogen_Mono) - (2 * H + O);
            vFragmentResults.push_back(fragment);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'T'] + ElementMassUtils::Hydrogen_Mono) - (3 * H + P + 4 * O);
            vFragmentResults.push_back(fragment);
        }
    }
    if (oligonucleotide.sSequence.find('C') != string::npos)
    {
        fragment.sSequence = "C";
        fragment.dMass = Param::g_staticParams.massUtility.pdBaseMass[(int)'C'] + ElementMassUtils::Hydrogen_Mono;
        if (Param::g_staticParams.options.bPositive)
            fragment.iCharge = 1;
        else 
            fragment.iCharge = -1;
      
        //M-base fragment
        for (int i = 1; i <= abs(iCharge); i++)
        {
            fragment.sSequence = oligonucleotide.sSequence;
            fragment.iCharge = i * iCharge / abs(iCharge);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'C'] + ElementMassUtils::Hydrogen_Mono);
            vFragmentResults.push_back(fragment);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'C'] + ElementMassUtils::Hydrogen_Mono) - (H + P + 3 * O);
            vFragmentResults.push_back(fragment);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'C'] + ElementMassUtils::Hydrogen_Mono) - (2 * H + O);
            vFragmentResults.push_back(fragment);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'C'] + ElementMassUtils::Hydrogen_Mono) - (3 * H + P + 4 * O);
            vFragmentResults.push_back(fragment);
        }
    }
    if (oligonucleotide.sSequence.find('G') != string::npos)
    {
        fragment.sSequence = "G";
        fragment.dMass = Param::g_staticParams.massUtility.pdBaseMass[(int)'G'] + ElementMassUtils::Hydrogen_Mono;
        if (Param::g_staticParams.options.bPositive)
            fragment.iCharge = 1;
        else 
            fragment.iCharge = -1;

        //M-base fragment
        for (int i = 1; i <= abs(iCharge); i++)
        {
            fragment.sSequence = oligonucleotide.sSequence;
            fragment.iCharge = i * iCharge / abs(iCharge);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'G'] + ElementMassUtils::Hydrogen_Mono);
            vFragmentResults.push_back(fragment);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'G'] + ElementMassUtils::Hydrogen_Mono) - (H + P + 3 * O);
            vFragmentResults.push_back(fragment);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'G'] + ElementMassUtils::Hydrogen_Mono) - (2 * H + O);
            vFragmentResults.push_back(fragment);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'G'] + ElementMassUtils::Hydrogen_Mono) - (3 * H + P + 4 * O);
            vFragmentResults.push_back(fragment);
        }
    }
    if (oligonucleotide.sSequence.find('U') != string::npos)
    {
        fragment.sSequence = "U";
        fragment.dMass = Param::g_staticParams.massUtility.pdBaseMass[(int)'U'] + ElementMassUtils::Hydrogen_Mono;
        if (Param::g_staticParams.options.bPositive)
            fragment.iCharge = 1;
        else 
            fragment.iCharge = -1;
      
        //M-base fragment
        for (int i = 1; i <= abs(iCharge); i++)
        {
            fragment.sSequence = oligonucleotide.sSequence;
            fragment.iCharge = i * iCharge / abs(iCharge);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'U'] + ElementMassUtils::Hydrogen_Mono);
            vFragmentResults.push_back(fragment);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'U'] + ElementMassUtils::Hydrogen_Mono) - (H + P + 3 * O);
            vFragmentResults.push_back(fragment);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'U'] + ElementMassUtils::Hydrogen_Mono) - (2 * H + O);
            vFragmentResults.push_back(fragment);

            fragment.dMass = oligonucleotide.dMass - (Param::g_staticParams.massUtility.pdBaseMass[(int)'U'] + ElementMassUtils::Hydrogen_Mono) - (3 * H + P + 4 * O);
            vFragmentResults.push_back(fragment);
        }
    }
    return vFragmentResults;
}

bool OligoProcessUtils::ReadNucleicAcidAndDigest(ThreadPool *tp, vector<NucleicAcid>& vNucleicAcid, vector<Oligonucleotide>& vOligonucleotide)
{
   ThreadPool *pSearchThreadPool = tp; //reuse the threadpool
   bool bSucceeded = true;

   DBInfo& databaseInfo = Param::g_staticParams.databaseInfo;
   bSucceeded = OligoProcessUtils::ReadNucleicAcid(databaseInfo, vNucleicAcid);
   if (!bSucceeded)
       return false;

   string sStaticModFilename = Param::g_staticParams.modificationInfo.sStaticModificationFilename;
   bSucceeded = ModificationUtils::AddStaticMod(vNucleicAcid, sStaticModFilename);
   if (!bSucceeded)
       return false;

   for (int index = 0; index < vNucleicAcid.size(); index++)
   {
       NucleicAcid& nucleicAcid = vNucleicAcid[index];
       while (pSearchThreadPool->jobs_.size() >= CommonValues::iMaxThreadProc)
       {
           pSearchThreadPool->wait_for_available_thread();
       }
       pSearchThreadPool->doJob(std::bind(DigestNucleicAcid, nucleicAcid, std::ref(vOligonucleotide)));
       
       if (!(index % 50))
       {
           WriteOutUtils::ShowProgressBar(index, vNucleicAcid.size(), "Digest Nucleic Acid", "");
       }

   }
   pSearchThreadPool->wait_on_threads();
   WriteOutUtils::ShowProgressBar(vNucleicAcid.size(), vNucleicAcid.size(), "Digest Nucleic Acid", "done\n");

    if (!bSucceeded)
        return false;
    else
        return true;
}

bool OligoProcessUtils::ReadNucleicAcid(DBInfo& databaseInfo, vector<NucleicAcid>& vNucleicAcid)
{
   bool bSucceeded = true;
   FILE *fp;
   int iTmpCh = 0;
   unsigned long long lEndPos = 0;
   unsigned long long lCurrPos = 0;
   bool bTrimDescr = false;
   
   databaseInfo.uliTotNTCount = 0;
   databaseInfo.iTotalNumNucleicAcid = 0;

   if ((fp=fopen(databaseInfo.sFilename.c_str(), "rb")) == NULL)
   {
        cout << " Error (1) - cannot read database file \"" + string(databaseInfo.sDatabase) + "\n";
        return false;
   }

   fseek(fp, 0, SEEK_END);
   lEndPos=ftell(fp);
   rewind(fp);

   // Load database entry header.
   lCurrPos = ftell(fp);
   iTmpCh = getc(fp);

   char szBuf[CommonValues::iFileBufferSize];
   int  iNumBadChars = 0; // count # of bad (non-printing) characters in header 
   bool bHeadOfFasta = true;
   // Loop through entire database.
   while(!feof(fp))
   {
      NucleicAcid nucleicAcid = NucleicAcid();
      if (bHeadOfFasta)
      {
         // skip through whitespace at head of line
         while (isspace(iTmpCh))
            iTmpCh = getc(fp);

         // skip comment lines
         if (iTmpCh == '#')
         {
            // skip to description line
            while ((iTmpCh != '\n') && (iTmpCh != '\r') && (iTmpCh != EOF))
               iTmpCh = getc(fp);
         }

         bHeadOfFasta = false;
      }

      if (iTmpCh == '>') // Expect a '>' for sequence header line.
      {
         // grab file pointer here for this sequence entry
         // this will be stored for protein references for each matched entry and
         // will be used to retrieve actual protein references when printing output
         nucleicAcid.lFilePosition = ftell(fp);

         bTrimDescr = false;
         while (((iTmpCh = getc(fp)) != '\n') && (iTmpCh != '\r') && (iTmpCh != EOF))
         {
            if (!bTrimDescr && iscntrl(iTmpCh))
               bTrimDescr = true;

            if (!bTrimDescr && nucleicAcid.sName.size() < (CommonValues::iMaxNameSize - 1))
            {
               if (iTmpCh < 32 || iTmpCh>126)  // sanity check for reading binary (index) file
               {
                  iNumBadChars++;
                  if (iNumBadChars > 20)
                  {
                     logerr(" Too many non-printing characters in database header lines; wrong file type/format?\n");
                     fclose(fp);
                     return false;
                  }
               }
               else
                  nucleicAcid.sName += iTmpCh;
            }
         }

         if (nucleicAcid.sName.length() <= 0)
         {
            string strErrorMsg = " Error - zero length sequence description; wrong database file/format?\n";
            // g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
            logerr(strErrorMsg.c_str());
            return false;
         }

         // Load sequence
         while (((iTmpCh=getc(fp)) != '>') && (iTmpCh != EOF))
         {
            if ('a' <= iTmpCh && iTmpCh <= 'z')
            {
               nucleicAcid.sSequence += iTmpCh - 32;  // convert toupper case so subtract 32 (i.e. 'A'-'a')
               databaseInfo.uliTotNTCount++;
            }
            else if ('A' <= iTmpCh && iTmpCh <= 'Z')
            {
               nucleicAcid.sSequence += iTmpCh;
               databaseInfo.uliTotNTCount++;
            }
            else if (iTmpCh == '*')  // stop codon
            {
               nucleicAcid.sSequence += iTmpCh;
            }
         }

         nucleicAcid.iWhichNucleicAcid = databaseInfo.iTotalNumNucleicAcid;
         databaseInfo.iTotalNumNucleicAcid++;
         vNucleicAcid.push_back(nucleicAcid);

         if (!(databaseInfo.iTotalNumNucleicAcid % 100))  // output procedure
         {
             lCurrPos = ftell(fp);
             WriteOutUtils::ShowProgressBar(lCurrPos, lEndPos, "     - Read Nucleic Acid progress: ", "");
         }

      }
      else
      {
         fgets(szBuf, sizeof(szBuf), fp);
         iTmpCh = getc(fp);
      }
   }

   {
       lCurrPos = ftell(fp);
       WriteOutUtils::ShowProgressBar(lCurrPos, lEndPos, "     - Read Nucleic Acid progress: ", "done\n");
   }

   fclose(fp);
   return true;
}

void OligoProcessUtils::DigestNucleicAcid(NucleicAcid nucleicAcid, vector<Oligonucleotide>& vOligonucleotide) //会使用多线程
{
   int iLenOligonucleotide = 0;
   int iLenNucleicAcid = nucleicAcid.sSequence.size();
   int iStartPos = 0;
   int iEndPos = 0;           // iEndPos == iStartPos at beginning

   // Search through entire protein.
   while (iStartPos < iLenNucleicAcid)
   {
      double dCalcMass = 0.0;

      iLenOligonucleotide = iEndPos - iStartPos + 1;  //include the iEndPos NT

      if (IntRange::inRange(Param::g_staticParams.options.oligonucleotideLengthRange, iLenOligonucleotide) // check length avaliable
         && CheckEnzymeTermini(nucleicAcid.sSequence, iStartPos, iEndPos, Param::g_staticParams.enzyme))   // check enzyme cut position
      {
         vector<Modification> vCandidate3end;
         vector<Modification> vCandidate5end;
         if (iStartPos == 0)
         {
            vCandidate5end.push_back(nucleicAcid.end5TermMod);
         }
         else
         {
            for (Modification & mod : Param::g_staticParams.modificationInfo.vAdd5terminusOligonucleotide)
            {
               vCandidate5end.push_back(mod);
            }
         }

         if (iEndPos == iLenNucleicAcid - 1)
         {  
            vCandidate3end.push_back(nucleicAcid.end3TermMod);
         }
         else
         {
            for (Modification & mod : Param::g_staticParams.modificationInfo.vAdd3terminusOligonucleotide)
            {
               vCandidate3end.push_back(mod);
            }
         }

         for (auto& _5endMod : vCandidate5end)
         {
            for (auto& _3endMod : vCandidate3end)
            {
               Oligonucleotide oligonucleotide;
               Oligonucleotide::NucleicAcidIndex whichNucleicAcid;
               whichNucleicAcid.iWhichNucleicAcid = nucleicAcid.iWhichNucleicAcid;
               whichNucleicAcid.iNABeginPos = iStartPos;
               whichNucleicAcid.iNAEndPos = iEndPos;
               oligonucleotide.sSequence = nucleicAcid.sSequence.substr(iStartPos, iLenOligonucleotide);
               oligonucleotide.whichNucleicAcid.push_back(whichNucleicAcid);
               oligonucleotide.end3TermMod = _3endMod;
               oligonucleotide.end5TermMod = _5endMod;

               for (const auto& mod : nucleicAcid.vMods) 
               {
                  int modPos = mod.first;
                  if (modPos >= iStartPos && modPos <= iEndPos) {
                     oligonucleotide.vMods.push_back({modPos - iStartPos, mod.second});
                  }
               }

               dCalcMass = MassUtils::CalculateOligonucleotideMass(oligonucleotide, true, true, -1); //no consider of various modification (just for original mass check)
               
               if (dCalcMass <= Param::g_staticParams.options.oligonucleotideMassRange.dEnd && oligonucleotide.sSequence.find('*') == std::string::npos)
               {
                   Threading::LockMutex(Param::g_vOligonucleotideMutex);
                   vOligonucleotide.push_back(oligonucleotide);
                   Threading::UnlockMutex(Param::g_vOligonucleotideMutex);
               }
            }
         }
      }

      // Increment end.
      if (dCalcMass <= Param::g_staticParams.options.oligonucleotideMassRange.dEnd && iEndPos < iLenNucleicAcid - 1 && iLenOligonucleotide < Param::g_staticParams.options.oligonucleotideLengthRange.iEnd)
         iEndPos++;
      // Increment start, reset end.
      else if (dCalcMass > Param::g_staticParams.options.oligonucleotideMassRange.dEnd || iEndPos == iLenNucleicAcid - 1 || iLenOligonucleotide == Param::g_staticParams.options.oligonucleotideLengthRange.iEnd)
      {
         iStartPos++;                           // Increment start of peptide.
         iEndPos = iStartPos;
      }
   }
   return;
}

bool OligoProcessUtils::CheckEnzymeTermini(string seq, int iStart, int iEnd, Enzyme& enzyme)
{
   if (enzyme.iEnzymeSelected)
   {
      bool bBeginCleavage = false;
      bool bEndCleavage = false;
      bool bBreakPoint;
      int iCountInternalCleavageSites=0;

      // todo: need fix according to the enzyme design (consider to form a class)
      bBeginCleavage = (iStart==0
            || seq[iStart-1] == '*'
            || ((enzyme.sEnzymeBreakNT.find(seq[iStart - enzyme.iSearchEnzymeOffSet]) != std::string::npos)
               && enzyme.sEnzymeNoBreakNT.find(seq[iStart - 1 - enzyme.iSearchEnzymeOffSet]) == std::string::npos));

      bEndCleavage = (iEnd==(int)(seq.size() - 1)
            || seq[iEnd +1 ] == '*'
            || ((enzyme.sEnzymeBreakNT.find(seq[iEnd + 1 - enzyme.iSearchEnzymeOffSet]) != std::string::npos)
               && enzyme.sEnzymeNoBreakNT.find(seq[iEnd - enzyme.iSearchEnzymeOffSet]) == std::string::npos));
            
      
      if (!(bBeginCleavage && bEndCleavage))
         return false;

      // Check number of missed cleavages count.
      // if N-term cleavage, iSearchEnzymeOffset==0, look at residues iStartPos+1 thru iEndPos to count missed cleavages
      // if C-term cleavage, iSearchEnzymeOffset==1, look at residues iStartPos thru iEndPos-1 to count missed cleavages
      int iBeginRef;
      int iEndRef;

      if (enzyme.iSearchEnzymeOffSet == 0)
      {
         iBeginRef = iStart + 1;
         iEndRef = iEnd;
      }
      else
      {
         iBeginRef = iStart;
         iEndRef = iEnd - 1;
      }

      for (int i = iBeginRef; i <= iEndRef; i++)
      {
         char cCurrentResidue = seq[i];
         char cFlankingResidue = (enzyme.iSearchEnzymeOffSet == 0 ? seq[i - 1] : seq[i + 1]);

         bBreakPoint = enzyme.sEnzymeBreakNT.find(cCurrentResidue) != std::string::npos
            && enzyme.sEnzymeNoBreakNT.find(cFlankingResidue) == std::string::npos;

         if (bBreakPoint)
         {
            if ((enzyme.iSearchEnzymeOffSet == 1 && i != iEnd)  // Ignore last residue.
                  || (enzyme.iSearchEnzymeOffSet == 0 && i != iStart))  // Ignore first residue.
            {
               iCountInternalCleavageSites++;
               if (iCountInternalCleavageSites > enzyme.iAllowedMissedCleavage)
                  return false;
            }
         }
      }
   }
   return true;
}

bool OligoProcessUtils::GenerateDecoyOligonucleotide(vector<Oligonucleotide>& vOligonucleotideList, int iDecoyOption)
{
    for (auto& oligonucleotide : vOligonucleotideList)
    {
        oligonucleotide.bDecoy = false;
        if (OligoProcessUtils::decoyDuplicateMap.find(oligonucleotide.sSequence) == OligoProcessUtils::decoyDuplicateMap.end())
        {
            OligoProcessUtils::decoyDuplicateMap[oligonucleotide.sSequence] = 1;
        }
    }

    vector<Oligonucleotide> vPreserveList;

    for (auto& oligonucleotide : vOligonucleotideList)
    {
        if (oligonucleotide.sSequence.size() <= 2)
            continue;
        Oligonucleotide decoy = oligonucleotide;
        decoy.bDecoy = true;
        decoy.iVarModificationIndex = 0;
        decoy.vMods.clear();

        for (int repeat_index = 0; repeat_index < 8; repeat_index++) // multi_shuffle to ensure decoy in enough
        {
            if (iDecoyOption == 0) //reverse
            {
                string reverseSequence = oligonucleotide.sSequence.substr(0, oligonucleotide.sSequence.size() - 1);
                reverse(reverseSequence.begin(), reverseSequence.end());
                decoy.sSequence = reverseSequence + oligonucleotide.sSequence.back();

                for (auto& mod : oligonucleotide.vMods)
                {
                    int index;
                    // 如果原本是最后一个字符的 modification，不进行位置反转
                    if (mod.first == oligonucleotide.sSequence.size() - 1)
                    {
                        index = mod.first;  // 保持最后一个字符的位置不变
                    }
                    else
                    {
                        // 对于其他字符，调整位置
                        index = oligonucleotide.sSequence.size() - 1 - mod.first;
                    }

                    // 将调整后的修改位置和修改内容加入到 decoy 的 vMods 中
                    decoy.vMods.push_back({ index, mod.second });
                }
            }
            else if (iDecoyOption == 1) //shuffle
            {
                int iSeqSize = oligonucleotide.sSequence.size();
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
                    reverseSequence.push_back(oligonucleotide.sSequence[index]);
                }
                reverseSequence.push_back(oligonucleotide.sSequence.back());
                decoy.sSequence = reverseSequence;

                // 更新 modifications 的位置
                for (const auto& mod : oligonucleotide.vMods)
                {
                    int originalPos = mod.first;
                    // 如果修改位置是最后一个字符，直接保持不变
                    if (originalPos == oligonucleotide.sSequence.size() - 1) {
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
            }

            if (OligoProcessUtils::decoyDuplicateMap.find(decoy.sSequence) == OligoProcessUtils::decoyDuplicateMap.end())
                vPreserveList.push_back(decoy);
        }
    }

    for (Oligonucleotide& oligonucleotide : vPreserveList)
    {
        vOligonucleotideList.push_back(oligonucleotide);
    }

    return true;

}
