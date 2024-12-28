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

#include "Param.h"
#include "MassUtils.h"

StaticParams                       Param::g_staticParams;         //static paramter, should not change in the program
MassChargeRange                    Param::g_mcRange;
vector<InputFile>                  Param::g_vInputFiles;
vector<NucleicAcid>                Param::g_vNucleicAcid;
vector<Oligonucleotide>            Param::g_vOligonucleotide;
vector<VarModification>            Param::g_vVarModification;
string                             Param::g_sCometVersion = comet_version;
//CometStatus                        Param::g_cometStatus;
Mutex                              Param::g_mcRangeMutex;
Mutex                              Param::g_vInputFileMutex;
Mutex                              Param::g_vNucleicAcidMutex;
Mutex                              Param::g_vOligonucleotideMutex;
Mutex                              Param::g_vVarModificationMutex;

//static function
static bool IsValidCometVersion(string sCometVersion)
{
    if (sCometVersion == "2024.0")
        return true;
    else
        return false;
}

static bool IsValidInputFile(string sFilename)
{
    ifstream file(sFilename);
    if (!file)
        return false;
    else
        return true;
}

static bool ValidateSequenceDatabaseFile()
{
    ifstream fsCheck;

    // if .idx database specified but does not exist, first see if corresponding
    // fasta exists and if it does, create the .idx file
    if (Param::g_staticParams.databaseInfo.sDatabase.find(".idx") != string::npos)
    {
        fsCheck = ifstream(Param::g_staticParams.databaseInfo.sDatabase);
        if (!fsCheck) //not .idx file
        {
            fsCheck = ifstream(Param::g_staticParams.databaseInfo.sFilename);
            if (!fsCheck)
            {
                string strErrorMsg = " Error - peptide index file \"" + Param::g_staticParams.databaseInfo.sDatabase + "\" and corresponding FASTA file\n         are both missing.\n";
                logerr(strErrorMsg.c_str());
                return false;
            }
            else
                return true;
        }
        else
        {
            fsCheck = ifstream(Param::g_staticParams.databaseInfo.sFilename);
            if (!fsCheck)
            {
                string strErrorMsg = " Error - peptide index file \"" + Param::g_staticParams.databaseInfo.sDatabase + "\" specified is present\n   but corresponding FASTA file\n         are both missing.\n";
                logerr(strErrorMsg.c_str());
                return false;
            }
            else
                return true;
        }
    }
    else
    {
        fsCheck = ifstream(Param::g_staticParams.databaseInfo.sDatabase);
        if (!fsCheck)
        {
            string strErrorMsg = " Error - Database file \"" + Param::g_staticParams.databaseInfo.sDatabase + "\" is missing.\n";
            logerr(strErrorMsg.c_str());
            return false;
        }
        else
        {
            return true;
        }
    }

    return true;
}

static bool ValidateScanRange()
{
    if (Param::g_staticParams.options.scanRange.iEnd < Param::g_staticParams.options.scanRange.iStart && Param::g_staticParams.options.scanRange.iEnd != 0)
    {
        string sErrorMsg = "Error - start scan is " + to_string(Param::g_staticParams.options.scanRange.iStart) + "but end scan is " + to_string(Param::g_staticParams.options.scanRange.iEnd) + "\n The end scan must be >= to the start scan.\n";
        logerr(sErrorMsg.c_str());
        return false;
    }

    return true;
}

static bool ValidatePeptideLengthRange()
{
    if (Param::g_staticParams.options.oligonucleotideLengthRange.iEnd < Param::g_staticParams.options.oligonucleotideLengthRange.iStart
        && Param::g_staticParams.options.oligonucleotideLengthRange.iEnd != 0)
    {
        char szErrorMsg[512];
        sprintf(szErrorMsg, " Error - peptide length range set as %d to %d.\n The maximum length must be >= to the minimum length.\n",
            Param::g_staticParams.options.oligonucleotideLengthRange.iStart,
            Param::g_staticParams.options.oligonucleotideLengthRange.iEnd);
        string strErrorMsg(szErrorMsg);
        // g_cometStatus.SetStatus(CometResult_Failed, strErrorMsg);
        logerr(szErrorMsg);
        return false;
    }

    return true;
}


// check tha validation of static parameters

bool Param::CheckParamsValidation()
{
    if (!ValidatePeptideLengthRange())
        return false;

    if (!ValidateScanRange())
        return false;

    if (!ValidateSequenceDatabaseFile)
        return false;

    return true;
}

void Param::initMutex()
{
    Threading::CreateMutex(&g_mcRangeMutex);

    Threading::CreateMutex(&g_vInputFileMutex);

    Threading::CreateMutex(&g_vNucleicAcidMutex);

    Threading::CreateMutex(&g_vOligonucleotideMutex);

    Threading::CreateMutex(&g_vVarModificationMutex);
}

void Param::destroyMutex()
{
    Threading::DestroyMutex(g_mcRangeMutex);

    Threading::DestroyMutex(g_vInputFileMutex);

    Threading::DestroyMutex(g_vNucleicAcidMutex);

    Threading::DestroyMutex(g_vOligonucleotideMutex);

    Threading::DestroyMutex(g_vVarModificationMutex);
}

bool Param::loadParam(string sParamFilename)
{
    string sLine;
    string sTmpStr;
    ifstream paramFile(sParamFilename);
    istringstream stringStream;

    if (!paramFile) // open fail
    {
        string sErrorMsg = "\n Comet version %s\n\n" + g_sCometVersion + " Error - cannot open parameter file \"%s\".\n" + sParamFilename;
        logerr(sErrorMsg.c_str());
        return false;
    }

    //get comet version
    getline(paramFile, sLine);
    stringStream = istringstream(sLine);
    stringStream >> sTmpStr >> sTmpStr >> sTmpStr;
    if (sTmpStr == "2024.02")
        g_sCometVersion = sTmpStr;
    else
    {
        string sErrorMsg = "\n Comet version %s\n\n" + g_sCometVersion + " The comet.params file is from version %s\n" + sTmpStr
            + " Please update your comet.params file.  You can generate\n" + " a new parameters file using \"comet -p\"\n\n";
        logerr(sErrorMsg.c_str());
        return false;
    }

    while (getline(paramFile, sLine)) {  // read parameter except enzyme
        if (sLine.find('#') != std::string::npos)
            sLine.erase(sLine.find('#'));

        if (sLine.find("[COMET_ENZYME_INFO]") != string::npos)
            break;

        stringStream = istringstream(sLine);
        string sParamName, sParamValue;
        getline(stringStream, sParamName, '=');
        getline(stringStream, sParamValue);

        sParamName.erase(std::remove_if(sParamName.begin(), sParamName.end(), [](char c) {return std::isspace(static_cast<char>(c));}), sParamName.end()); //去除两端多余的空格

        sParamValue.erase(std::remove_if(sParamValue.begin(), sParamValue.end(), [](char c) {return std::isspace(static_cast<char>(c));}), sParamValue.end());

        if (sParamName == "")
            continue;
        else if (sParamName == "database_name")
        {
            if (sParamValue.find(".idx") != string::npos)
            {
                g_staticParams.databaseInfo.sDatabase = sParamValue;
                g_staticParams.databaseInfo.sFilename = sParamValue.substr(0, sParamValue.find(".idx"));
            }
            else
            {
                g_staticParams.databaseInfo.sDatabase = sParamValue + ".idx";
                g_staticParams.databaseInfo.sFilename = sParamValue;
            }
        }
        else if (sParamName == "decoy_search")
        {
            g_staticParams.options.iDecoySearch = stoi(sParamValue);
        }
        else if (sParamName == "decoy_kind")
        {
            g_staticParams.options.iDecoyKind = stoi(sParamValue);
        }
        else if (sParamName == "minimum_xcorr")
        {
            g_staticParams.options.dMinimumXcorr = stod(sParamValue);
        }
//charge
        else if (sParamName == "max_fragment_charge")
        {
            g_staticParams.options.fragmentChargeRange.iEnd = stoi(sParamValue);
        }
        else if (sParamName == "max_precursor_charge")
        {
            g_staticParams.options.precursorChargeRange.iEnd = stoi(sParamValue);
        }
        else if (sParamName == "min_fragment_charge")
        {
            g_staticParams.options.fragmentChargeRange.iStart = stoi(sParamValue);
        }
        else if (sParamName == "min_precursor_charge")
        {
            g_staticParams.options.precursorChargeRange.iStart = stoi(sParamValue);
        }
        else if (sParamName == "positive_charge")
        {
            if (stoi(sParamValue))
                g_staticParams.options.bPositive = true;
            else
                g_staticParams.options.bPositive = false;
        }
        else if (sParamName == "adduct_ion")
        {
            g_staticParams.options.sAdductIon = sParamValue;
        }
//thread
        else if (sParamName == "num_threads")
        {
            g_staticParams.options.iNumThreads = stoi(sParamValue);
        }
//fragment ion
        else if (sParamName == "theoretical_fragment_ions")
        {
            g_staticParams.ionInformation.iTheoreticalFragmentIons = stoi(sParamValue);
        }
        else if (sParamName == "use_a_ions")
        {
            g_staticParams.ionInformation.vSelectedIonSeries.push_back(IonSeries::ion_a);
        }
        else if (sParamName == "use_b_ions")
        {
            g_staticParams.ionInformation.vSelectedIonSeries.push_back(IonSeries::ion_b);
        }
        else if (sParamName == "use_c_ions")
        {
            g_staticParams.ionInformation.vSelectedIonSeries.push_back(IonSeries::ion_c);
        }
        else if (sParamName == "use_d_ions")
        {
            g_staticParams.ionInformation.vSelectedIonSeries.push_back(IonSeries::ion_d);
        }
        else if (sParamName == "use_w_ions")
        {
            g_staticParams.ionInformation.vSelectedIonSeries.push_back(IonSeries::ion_w);
        }
        else if (sParamName == "use_x_ions")
        {
            g_staticParams.ionInformation.vSelectedIonSeries.push_back(IonSeries::ion_x);
        }
        else if (sParamName == "use_y_ions")
        {
            g_staticParams.ionInformation.vSelectedIonSeries.push_back(IonSeries::ion_y);
        }
        else if (sParamName == "use_z_ions")
        {
            g_staticParams.ionInformation.vSelectedIonSeries.push_back(IonSeries::ion_z);
        }
        else if (sParamName == "use_a_B_ions")
        {
            g_staticParams.ionInformation.vSelectedIonSeries.push_back(IonSeries::ion_a_B);
        }
        else if (sParamName == "use_NL_ions")
        {
            g_staticParams.ionInformation.vSelectedIonSeries.push_back(IonSeries::ion_NL);
        }
        else if (sParamName == "fragment_bin_tol")
        {
            g_staticParams.dInverseBinWidth = 1.0 / stod(sParamValue);
        }
        else if (sParamName == "fragment_bin_offset")
        {
            g_staticParams.dOneMinusBinOffset = 1.0 - stod(sParamValue);
        }
//tolerance
        else if (sParamName == "precursor_tolerance_upper")
        {
            g_staticParams.tolerances.precursorTolerance.dEnd = stod(sParamValue);
        }
        else if (sParamName == "precursor_tolerance_lower")
        {
            g_staticParams.tolerances.precursorTolerance.dStart = stod(sParamValue);
        }
        else if (sParamName == "precursor_tolerance_type")
        {
            g_staticParams.tolerances.iPrecursorToleranceType = stoi(sParamValue);
        }
        else if (sParamName == "precursor_tolerance_units")
        {
            g_staticParams.tolerances.iPrecursorToleranceUnits = stoi(sParamValue);
        }
        else if (sParamName == "fragment_tolerance_upper")
        {
            g_staticParams.tolerances.fragmentTolerance.dEnd = stod(sParamValue);
        }
        else if (sParamName == "fragment_tolerance_lower")
        {
            g_staticParams.tolerances.fragmentTolerance.dStart = stod(sParamValue);
        }
        else if (sParamName == "fragment_tolerance_type")
        {
            g_staticParams.tolerances.iFragmentToleranceType = stoi(sParamValue);
        }
        else if (sParamName == "fragment_tolerance_units")
        {
            g_staticParams.tolerances.iFragmentToleranceUnits = stoi(sParamValue);
        }
        else if (sParamName == "isotope_error")
        {
            g_staticParams.tolerances.iIsotopeSearchType = stoi(sParamValue);
        }
        else if (sParamName == "mass_type")
        {
            g_staticParams.massUtility.bMonoMass = stoi(sParamValue);
        }
        else if (sParamName == "deoxidation")
        {
            if (stoi(sParamValue))
                g_staticParams.massUtility.bDeoxidation = true;
            else
                g_staticParams.massUtility.bDeoxidation = false;
        }
//spectrum preprocess
        else if (sParamName == "remove_precursor_peak")
        {
        }
        else if (sParamName == "remove_precursor_tolerance")
        {
        }
        else if (sParamName == "clear_mz_range")
        {
            DoubleRange range;
            if (sParamValue.find('-') != string::npos)
            {
                range.dStart = stod(sParamValue.substr(0, sParamValue.find('-')));
                range.dEnd = stod(sParamValue.substr(sParamValue.find('-') + 1));
            }
            g_staticParams.options.clearMzRange = range;
        }
        else if (sParamName == "percentage_base_peak")
        {
            g_staticParams.options.dMinPercentageIntensity = stod(sParamValue);
        }
        else if (sParamName == "scan_range")
        {
            IntRange range;
            if (sParamValue.find('-') != string::npos)
            {
                range.iStart = stoi(sParamValue.substr(0, sParamValue.find('-')));
                range.iEnd = stoi(sParamValue.substr(sParamValue.find('-') + 1));
            }
            g_staticParams.options.scanRange = range;
        }
        else if (sParamName == "spectrum_batch_size")
        {
            g_staticParams.options.iSpectrumBatchSize = stoi(sParamValue);
        }
        else if (sParamName == "minimum_peaks")
        {
            g_staticParams.options.iMinPeaks = stoi(sParamValue);
        }
        else if (sParamName == "ms_level")
        {
            g_staticParams.options.iMSLevel = stoi(sParamValue);
        }
        else if (sParamName == "activation_method")
        {
        }
        else if (sParamName == "minimum_intensity")
        {
            g_staticParams.options.dMinIntensity = stod(sParamValue);
        }
//output
//modification
        else if (sParamName == "add_3term_oligonucleotide")
        {

            if (sParamValue.find("OH") != string::npos)
            {
                Modification mod;
                mod.sModName = "OH";
                mod.dModMass = ElementMassUtils::Oxygen_Mono + ElementMassUtils::Hydrogen_Mono; //理论上这里应该判断是否mono，但是懒
                g_staticParams.modificationInfo.vAdd3terminusOligonucleotide.push_back(mod);
            }

            if (sParamValue.find("cP") != string::npos)
            {
                Modification mod;
                mod.sModName = "cP";
                mod.dModMass = -1 * (ElementMassUtils::Phosphorus_Mono + 2 * ElementMassUtils::Oxygen_Mono);
                g_staticParams.modificationInfo.vAdd3terminusOligonucleotide.push_back(mod);
            }

            if (g_staticParams.modificationInfo.vAdd3terminusOligonucleotide.size() == 0)
            {
                Modification mod;
                mod.sModName = "";
                mod.dModMass = 0.0;
                g_staticParams.modificationInfo.vAdd3terminusOligonucleotide.push_back(mod);
            }

        }
        else if (sParamName == "add_5term_oligonucleotide")
        {
            if (sParamValue.find("OH") != string::npos)
            {
                Modification mod;
                mod.sModName = "OH";
                mod.dModMass = ElementMassUtils::Hydrogen_Mono;                                 //理论上这里应该判断是否mono，但是懒
                g_staticParams.modificationInfo.vAdd5terminusOligonucleotide.push_back(mod);
            }

            if (sParamValue.find("P") != string::npos)
            {
                Modification mod;
                mod.sModName = "P";
                mod.dModMass = ElementMassUtils::Phosphorus_Mono + 3 * ElementMassUtils::Oxygen_Mono + 2 * ElementMassUtils::Hydrogen_Mono;
                g_staticParams.modificationInfo.vAdd5terminusOligonucleotide.push_back(mod);
            }

            if (g_staticParams.modificationInfo.vAdd5terminusOligonucleotide.size() == 0)
            {
                Modification mod;
                mod.sModName = "";
                mod.dModMass = 0.0;
                g_staticParams.modificationInfo.vAdd5terminusOligonucleotide.push_back(mod);
            }
        }
        else if (sParamName == "add_3term_nucleic_acid")
        {
            if (sParamValue.find("OH") != string::npos)
            {
                Modification mod;
                mod.sModName = "OH";
                mod.dModMass = ElementMassUtils::Oxygen_Mono + ElementMassUtils::Hydrogen_Mono; //理论上这里应该判断是否mono，但是懒
                g_staticParams.modificationInfo.vAdd3terminusNucleicAcid.push_back(mod);
            }

            if (sParamValue.find("cP") != string::npos)
            {
                Modification mod;
                mod.sModName = "cP";
                mod.dModMass = -1 * (ElementMassUtils::Phosphorus_Mono + 2 * ElementMassUtils::Oxygen_Mono);
                g_staticParams.modificationInfo.vAdd3terminusNucleicAcid.push_back(mod);
            }

            if (g_staticParams.modificationInfo.vAdd3terminusNucleicAcid.size() == 0)
            {
                Modification mod;
                mod.sModName = "";
                mod.dModMass = 0.0;
                g_staticParams.modificationInfo.vAdd3terminusNucleicAcid.push_back(mod);
            }
        }
        else if (sParamName == "add_5term_nucleic_acid")
        {
            if (sParamValue.find("OH") != string::npos)
            {
                Modification mod;
                mod.sModName = "OH";
                mod.dModMass = ElementMassUtils::Hydrogen_Mono;                                 //理论上这里应该判断是否mono，但是懒
                g_staticParams.modificationInfo.vAdd5terminusNucleicAcid.push_back(mod);
            }

            if (sParamValue.find("P") != string::npos)
            {
                Modification mod;
                mod.sModName = "P";
                mod.dModMass = ElementMassUtils::Phosphorus_Mono + 3 * ElementMassUtils::Oxygen_Mono + 2 * ElementMassUtils::Hydrogen_Mono;
                g_staticParams.modificationInfo.vAdd5terminusNucleicAcid.push_back(mod);
            }

            if (g_staticParams.modificationInfo.vAdd5terminusNucleicAcid.size() == 0)
            {
                Modification mod;
                mod.sModName = "";
                mod.dModMass = 0.0;
                g_staticParams.modificationInfo.vAdd5terminusNucleicAcid.push_back(mod);
            }
        }
        else if (sParamName == "max_variable_mods_in_oligonucleotide")
        {
            g_staticParams.options.iMaxVarModPerOligonucleoitde = stoi(sParamValue);
        }
        else if (sParamName == "variable_modslist_file")
        {
            g_staticParams.modificationInfo.sVariableModificationFilename = sParamValue;
        }
        else if (sParamName == "static_modslist_file")
        {
            g_staticParams.modificationInfo.sStaticModificationFilename = sParamValue;
        }
//enzyme
        else if (sParamName == "search_enzyme_number")
        {
            g_staticParams.enzyme.iEnzymeSelected = stoi(sParamValue);
        }
        else if (sParamName == "allowed_missed_cleavage")
        {
            g_staticParams.enzyme.iAllowedMissedCleavage = stoi(sParamValue);
        }
// Oligonucleotide && Nucleic Acid Param
        else if (sParamName == "oligonucleotide_length_range")
        {
            IntRange range;
            if (sParamValue.find('-') != string::npos)
            {
                range.iStart = stoi(sParamValue.substr(0, sParamValue.find('-')));
                range.iEnd = stoi(sParamValue.substr(sParamValue.find('-') + 1));
            }
            g_staticParams.options.oligonucleotideLengthRange = range;
        }
        else if (sParamName == "digest_mass_range")
        {
            DoubleRange range;
            if (sParamValue.find('-') != string::npos)
            {
                range.dStart = stod(sParamValue.substr(0, sParamValue.find('-')));
                range.dEnd = stod(sParamValue.substr(sParamValue.find('-') + 1));
            }
            g_staticParams.options.oligonucleotideMassRange = range;
        }
        else
        {
            string sErrorMsg = " Warning - invalid parameter found: %s.  Parameter will be ignored.\n" + sParamName;
            logout(sErrorMsg.c_str());
        }

    }

    while (getline(paramFile, sLine))
    {
        stringStream = istringstream(sLine);
        stringStream >> sTmpStr;
        int enzymeNumber = stoi(sTmpStr.substr(0, sTmpStr.size() - 1));
        if (enzymeNumber == g_staticParams.enzyme.iEnzymeSelected)
        {
            stringStream >> g_staticParams.enzyme.sEnzymeName;
            stringStream >> sTmpStr;
            g_staticParams.enzyme.iSearchEnzymeOffSet = stoi(sTmpStr);
            stringStream >> g_staticParams.enzyme.sEnzymeBreakNT;
            stringStream >> g_staticParams.enzyme.sEnzymeNoBreakNT;
        }
    }

    MassUtils::AssignMass(Param::g_staticParams.massUtility);

    return true;
}

void Param::printParams(string sParamFilename, int option)
{
    ofstream paramFile(sParamFilename);
    if (!paramFile)
    {
        string sErrorMsg = "\n Comet version %s\n\n" + g_sCometVersion + " Error - cannot write file comet.params.new\n" + sParamFilename;
        logerr(sErrorMsg.c_str());
        exit(0);
    }

    // 1. 写入基本的注释和 Comet 版本信息
    paramFile << "# comet_version " + g_sCometVersion << endl;
    paramFile << "# Comet MS/MS search engine parameters file." << endl;
    paramFile << "# Everything following the '#' symbol is treated as a comment." << endl;

    // 2. 数据库相关设置
    paramFile << "#" << endl;
    paramFile << "# database parameters" << endl;
    paramFile << "#" << endl;
    paramFile << "database_name = xxx.fasta.idx" << endl;
    paramFile << "decoy_search = 1                                      # 0=no (default), 1=decoy" << endl;
    paramFile << "decoy_kind = 1                                        # 0=reverse, 1=shuffle" << endl;

    // 3. 质量和容差相关设置
    paramFile << "#" << endl;
    paramFile << "# masses && tolerance" << endl;
    paramFile << "#" << endl;
    paramFile << "precursor_tolerance_upper = 16.0                    # upper bound of the precursor mass tolerance" << endl;
    paramFile << "precursor_tolerance_lower = -16.0                   # lower bound of the precursor mass tolerance; USUALLY NEGATIVE TO BE LOWER THAN 0" << endl;
    paramFile << "precursor_tolerance_units = 2                       # 0=amu, 1=mmu, 2=ppm" << endl;
    paramFile << "precursor_tolerance_type = 1                        # 0=MH+ (default), 1=precursor m/z; only valid for amu/mmu tolerances" << endl;
    paramFile << "fragment_tolerance_upper = 16.0                     # upper bound of the precursor mass tolerance" << endl;
    paramFile << "fragment_tolerance_lower = -16.0                    # lower bound of the precursor mass tolerance; USUALLY NEGATIVE TO BE LOWER THAN 0" << endl;
    paramFile << "fragment_tolerance_units = 2                        # 0=amu, 1=mmu, 2=ppm" << endl;
    paramFile << "fragment_tolerance_type = 1                         # 0=MH+ (default), 1 = m/z; only valid for amu/mmu tolerances" << endl;
    paramFile << "isotope_error = 0                                   # 0=off, 1=0/1 (C13 error), 2=0/1/2, 3=0/1/2/3, 4=-1/0/1/2/3, 5=-1/0/1" << endl;
    paramFile << "mass_type = 1                                       # 0=average masses, 1=monoisotopic masses" << endl;
    paramFile << "deoxidation = 0                                     # 0 = RNA, 1 = DNA" << endl;

    // 4. 电荷状态相关设置
    paramFile << "#" << endl;
    paramFile << "# charge" << endl;
    paramFile << "#" << endl;
    paramFile << "max_fragment_charge = 3                             # set maximum fragment charge state to analyze" << endl;
    paramFile << "max_precursor_charge = 6                            # set maximum precursor charge state to analyze" << endl;
    paramFile << "min_fragment_charge = 1" << endl;
    paramFile << "min_precursor_charge = 1" << endl;
    paramFile << "positive_charge = 1                                 # 0 = negtive, 1 = positive" << endl;

    // 5. 线程设置
    paramFile << "#" << endl;
    paramFile << "# thread parameter" << endl;
    paramFile << "#" << endl;
    paramFile << "num_threads = 0                                   # 0=poll CPU to set num threads; else specify num threads directly (recommand max 128)" << endl;

    // 6. 酶切设定
    paramFile << "#" << endl;
    paramFile << "# search enzyme parameter" << endl;
    paramFile << "#" << endl;
    paramFile << "search_enzyme_number = 2                # choose from list at end of this params file" << endl;
    paramFile << "allowed_missed_cleavage = 2             # maximum value is 5; for enzyme search" << endl;

    // 7. 各种参数和静态修改
    paramFile << "#" << endl;
    paramFile << "# various && static modifications" << endl;
    paramFile << "#" << endl;
    paramFile << "add_3term_oligonucleotide = OH                                              # option = [OH, cP], can choose multiple" << endl;
    paramFile << "add_5term_oligonucleotide = OH                                              # option = [OH, P]" << endl;
    paramFile << "add_3term_nucleic_acid = OH                                                 # option = [OH, cP]" << endl;
    paramFile << "add_5term_nucleic_acid = OH                                                 # option = [OH, P]" << endl;
    paramFile << "max_variable_mods_in_oligonucleotide = 5                                    # limit variable mods to subset of specified proteins if this file is specified & present" << endl;
    paramFile << "variable_modslist_file = VarMod.params" << endl;
    paramFile << "static_modslist_file = StaticMod.params" << endl;

    // 8. 碎片离子参数
    paramFile << "#" << endl;
    paramFile << "# fragment ions" << endl;
    paramFile << "#" << endl;
    paramFile << "fragment_bin_tol = 0.02                # binning to use on fragment ions" << endl;
    paramFile << "fragment_bin_offset = 0.0              # offset position to start the binning (0.0 to 1.0)" << endl;
    paramFile << "use_a_ions = 0" << endl;
    paramFile << "use_b_ions = 1" << endl;
    paramFile << "use_c_ions = 0" << endl;
    paramFile << "use_d_ions = 0" << endl;
    paramFile << "use_w_ions = 0" << endl;
    paramFile << "use_x_ions = 0" << endl;
    paramFile << "use_y_ions = 1" << endl;
    paramFile << "use_z_ions = 0" << endl;
    paramFile << "use_a_B_ions = 0" << endl;
    paramFile << "use_NL_ions = 0                        # 0=no, 1=yes to consider NH3/H2O neutral loss peaks" << endl;

    // 9. Oligonucleotide 和核酸参数
    paramFile << "#" << endl;
    paramFile << "# Oligonucleotide && Nucleic Acid Param" << endl;
    paramFile << "#" << endl;
    paramFile << "digest_mass_range = 200.0-5000.0                # MH+ peptide mass range to analyze" << endl;
    paramFile << "fragment_mass_range = 200.0-5000.0              # fragment mass range to analyze" << endl;
    paramFile << "oligonucleotide_length_range = 3-50             # minimum and maximum peptide length to analyze (default min 1 to allowed max 51)" << endl;

    // 10. output 参数
    paramFile << "#" << endl;
    paramFile << "# output parameters" << endl;
    paramFile << "#" << endl;
    paramFile << "num_results = 100                               # number of results to store internally for Sp rank only; if Sp rank is not used, set this to num_output_lines" << endl;

    // 10. 光谱处理参数
    paramFile << "#" << endl;
    paramFile << "# spectral processing" << endl;
    paramFile << "#" << endl;
    paramFile << "scan_range = 0-0                       # start and end scan range to search; either entry can be set independently" << endl;
    paramFile << "minimum_peaks = 10                     # required minimum number of peaks in spectrum to search (default 10)" << endl;
    paramFile << "minimum_intensity = 0                 # minimum intensity value to read in" << endl;
    paramFile << "remove_precursor_peak = 0              # 0=no, 1=yes, 2=all charge reduced precursor peaks (for ETD), 3=phosphate neutral loss peaks" << endl;
    paramFile << "remove_precursor_tolerance = 1.5       # +- Da tolerance for precursor removal" << endl;
    paramFile << "clear_mz_range = 0.0-0.0               # clear out all peaks in the specified m/z range e.g. remove reporter ion region of TMT spectra" << endl;
    paramFile << "percentage_base_peak = 0.0             # specify a percentage (e.g. \"0.05\" for 5%) of the base peak intensity as a minimum intensity threshold" << endl;
    paramFile << "ms_level = 2                           # MS level to analyze, valid are levels 2 (default) or 3" << endl;
    paramFile << "activation_method = ALL                # activation method; used if activation method set; allowed ALL, CID, ECD, ETD, ETD+SA, PQD, HCD, IRMPD, SID" << endl;
    paramFile << "spectrum_batch_size = 15000            # max. # of spectra to search at a time; 0 to search the entire scan range in one loop" << endl;

    // 11. 酶信息（必须放在文件末尾）
    paramFile << "#" << endl;
    paramFile << "# COMET_ENZYME_INFO _must_ be at the end of this parameters file" << endl;
    paramFile << "# Enzyme entries can be added/deleted/edited" << endl;
    paramFile << "# format: <enzyme name> <whether 3' cut> <cut nucleotide> <no cut nucleotide>" << endl;
    paramFile << "#" << endl;
    paramFile << "[COMET_ENZYME_INFO]" << endl;
    paramFile << "0. Cut_everywhere 0 - -" << endl;
    paramFile << "1. Enzyme_A 0 A -" << endl;
    paramFile << "2. Enzyme_T1 1 G  -" << endl;
    paramFile << "3. No_cut 1 @ @" << endl;

    logout("\n Created:  comet.params.new\n\n");

    exit(1);
}

bool Param::ProcessCmdLine(int argc, char* argv[])
{
    vector<string> vArgv;
    int iPrintParams = 0;
    int iArgvIndex = 1;
    string sTmpStr;
    string sParamsFilename;

    if (argc == 1)
    {
        string sErrorMsg = "\n Comet version %s\n\n" + g_sCometVersion + " Error - no input files specified so nothing to do.\n";
        logerr(sErrorMsg.c_str());
        exit(1);
    }

    sParamsFilename = "./comet.params";                    //default parameter file name, may change in SetOptions(...) function

    //read all args
    while (iArgvIndex < argc)
    {
        sTmpStr = string(argv[iArgvIndex]);
        vArgv.push_back(sTmpStr);
        ++iArgvIndex;
    }

    //reload parameter file
    for (string& arg : vArgv)
    {
        if (arg[0] == '-' && arg[1] == 'P')
            sParamsFilename = arg.substr(2);
        if (arg[0] == '-' && arg[1] == 'p')
            Param::printParams("./comet.params.new", 1);
        if (arg[0] == '-' && arg[1] == 'q')
            Param::printParams("./comet.params.new", 2);
    }

    if (!Param::loadParam(sParamsFilename))
    {
        string sErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n" + " Error - load parameter file fail. \n";
        logerr(sErrorMsg.c_str());
        return false;
    }

    //set options and read inputfile
    for (string& arg : vArgv)
    {
        if (arg[0] == '-')
            Param::SetOptions(arg);
        else
        {
            if (!IsValidInputFile(arg))
            {
                string sErrorMsg = "\n Comet version " + g_sCometVersion + "\n\n" + " Error - input file \"" + arg + "\" not found.\n";
                logerr(sErrorMsg.c_str());
                return false;
            }
            InputFile inputFile;
            inputFile.sFileName = arg;
            inputFile.iFirstScan = Param::g_staticParams.options.scanRange.iStart;
            inputFile.iLastScan = Param::g_staticParams.options.scanRange.iEnd;
            Param::g_vInputFiles.push_back(inputFile);
        }
    }

    return true;
}

bool Param::printUsage(string sCmd)
{
    string sLogMsg;

    logout("\n");
    sLogMsg = " Comet version \" " + g_sCometVersion + "\"\n";
    logout(sLogMsg.c_str());
    logout("\n");
    sLogMsg = " Comet usage:  " + sCmd + "\"\n";
    logout(sLogMsg.c_str());
    logout("\n");
    logout(" Supported input formats include mzXML, mzML, Thermo raw, mgf, and ms2 variants (cms2, bms2, ms2)\n");

    logout("\n");
    logout("       options:  -p         to print out a comet.params.new file\n");
    logout("                 -q         to print out a comet.params.new file with more parameter entries\n");
    logout("                 -P<params> to specify an alternate parameters file (default comet.params)\n");
    logout("                 -D<dbase>  to specify a sequence database, overriding entry in parameters file\n");
    logout("                 -F<num>    to specify the first/start scan to search, overriding entry in parameters file\n");
    logout("                 -L<num>    to specify the last/end scan to search, overriding entry in parameters file\n");
    logout("                            (-L option is required if -F option is used)\n");
    logout("\n");
    sLogMsg = "       example:  " + sCmd + "file1.mzXML file2.mzXML\n";
    logout(sLogMsg.c_str());
    sLogMsg = "            or   " + sCmd + "-F1000 -L1500 file1.mzXML    <- to search scans 1000 through 1500\n";
    logout(sLogMsg.c_str());
    sLogMsg = "            or   " + sCmd + "-PParams.txt *.mzXML         <- use parameters in the file 'Params.txt'\n";
    logout(sLogMsg.c_str());
    logout("\n");

    exit(1);
}

bool Param::SetOptions(string arg)
{
    switch (arg[1])
    {
    case 'D':                        // Alternate sequence database.
        if (arg.size() <= 2)
            logerr("Missing text for parameter option -D<database>.  Ignored.\n");
        else
            Param::g_staticParams.databaseInfo.sFilename = arg.substr(2);
        break;
    case 'F':   // first scan
        if (arg.size() <= 2)
            logerr("Missing text for parameter option -F<num>.  Ignored.\n");
        else
            Param::g_staticParams.options.scanRange.iStart = stoi(arg.substr(2));
        break;
    case 'L':  // last scan
        if (arg.size() <= 2)
            logerr("Missing text for parameter option -L<num>.  Ignored.\n");
        else
            Param::g_staticParams.options.scanRange.iEnd = stoi(arg.substr(2));
        break;
    case 'B':    // batch size
        if (arg.size() <= 2)
            logerr("Missing text for parameter option -B<num>.  Ignored.\n");
        else
            Param::g_staticParams.options.iSpectrumBatchSize = stoi(arg.substr(2));
        break;
    default:
        break;
    }
    return true;
}













