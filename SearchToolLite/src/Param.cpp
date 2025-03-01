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

#include <string>
#include <map>
#include <fstream>
#include <sstream>

// ConfigParser
bool ConfigParser::parse(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false; // Return false if the file cannot be opened
    }
    std::string line;
    std::string current_section;

    while (std::getline(file, line)) {
        // Remove leading and trailing whitespace
        line.erase(0, line.find_first_not_of(" \t"));
        line.erase(line.find_last_not_of(" \t") + 1);

        // Remove comments
        size_t comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }

        // Skip empty lines
        if (line.empty()) {
            continue;
        }

        // Check if this is a section header
        if (line[0] == '[' && line[line.length() - 1] == ']') {
            current_section = line.substr(1, line.length() - 2);
            continue;
        }

        // Parse key-value pair
        size_t delimiter_pos = line.find('=');
        if (delimiter_pos != std::string::npos) {
            std::string key = line.substr(0, delimiter_pos);
            std::string value = line.substr(delimiter_pos + 1);

            // Remove whitespace
            key.erase(0, key.find_first_not_of(" \t"));
            key.erase(key.find_last_not_of(" \t") + 1);
            value.erase(0, value.find_first_not_of(" \t"));
            value.erase(value.find_last_not_of(" \t") + 1);

            config[current_section][key] = value;
        }
    }
    return true;
}

std::string ConfigParser::get(const std::string& section, const std::string& option) {
    if (config.find(section) != config.end()) {
        if (config[section].find(option) != config[section].end()) {
            return config[section][option];
        }
        else
            throw std::runtime_error("Option " + option + " not found in section " + section + " in config file");
    }
    else
        throw std::runtime_error("Section " + section + " not found in config file");
}

// TSVParser
void TSVParser::parseHeaders(const std::string& headerLine) {
    std::stringstream ss(headerLine);
    std::string header;
    while (std::getline(ss, header, '\t')) {
        header.erase(0, header.find_first_not_of(" \t"));
        header.erase(header.find_last_not_of(" \t") + 1);
        headers.push_back(header);
    }
}

std::map<std::string, std::string> TSVParser::parseLine(const std::string& line) {
    std::map<std::string, std::string> row;
    std::stringstream ss(line);
    std::string value;
    size_t columnIndex = 0;

    while (std::getline(ss, value, '\t')) {
        value.erase(0, value.find_first_not_of(" \t"));
        value.erase(value.find_last_not_of(" \t") + 1);
        if (columnIndex < headers.size()) {
            row[headers[columnIndex]] = value;
        }
        columnIndex++;
    }

    return row;
}

bool TSVParser::readFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }
    bool bFirstLine = true;
    std::string line;
    // Read and parse data lines
    while (std::getline(file, line)) {
        // Remove comments
        size_t comment_pos = line.find('#');
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }
        // Skip empty lines
        if (line.empty()) {
            continue;
        }
        // Parse headers
        if (bFirstLine) {
            parseHeaders(line);
            bFirstLine = false;
        }
        // Parse data
        else {
            data.push_back(parseLine(line));
        }
    }

    return true;
}

std::vector<std::map<std::string, std::string>> TSVParser::getData() const {
    return data;
}

std::vector<std::string> TSVParser::getHeaders() const {
    return headers;
}

//Options

Options::Options()
{
    sDatabaseName = "./resource/ecoli.fasta";
    bDecoySearch = true; // 1 = true
    iDecoyKind = 1; // 1 = shuffle

    precursorTolerance = DoubleRange(-16.0, 16.0);
    fragmentTolerance = DoubleRange(-16.0, 16.0);
    iIsotopeError = 0;
    bDeoxidation = false; // 0 = false

    fragmentChargeRange = IntRange(1, 3);
    precursorChargeRange = IntRange(1, 6);
    bPositiveCharge = false; // 0 = false
    sAdductIon = "@H+";

    iThreadNum = 2;

    sEnzymeDatabase = "./resource/enzyme.ini";
    sEnzymeName = "enzyme_A";
    iAllowedMissed = 2;

    sVarModsDatabase = "./resource/varMod.ini";
    sStaticModsDatabase = "./resource/staticMod.ini";
    sOligo3Mod = "OH,cP";
    sOligo5Mod = "OH,P";
    sNucleic3Mod = "OH,cP";
    sNucleic5Mod = "OH,P";
    iMaxVarMods = 5;

    bUseAIons = true; // 1 = true
    bUseBIons = true; // 1 = true
    bUseCIons = true; // 1 = true
    bUseDIons = true; // 1 = true
    bUseWIons = true; // 1 = true
    bUseXIons = true; // 1 = true
    bUseYIons = true; // 1 = true
    bUseZIons = true; // 1 = true
    bUseABIons = true; // 1 = true
    bUseNLIons = true; // 1 = true

    digestMassRange = DoubleRange(200.0, 5000.0);
    fragmentMassRange = DoubleRange(200.0, 5000.0);
    oligoLengthRange = IntRange(3, 50);

    scanRange = IntRange(0, 0);
    dMinimumIntensity = 0;
    bRemovePrecursorPeak = false; // 0 = false
    dRemovePrecursorTolerance = 1.5;
    clearMZRange = DoubleRange(0.0, 0.0);
    dPercentageBasePeak = 0.0;
    iMSLevel = 2;
    iSpectrumBatchSize = 15000;
}

bool Options::loadParameters(const std::string& filename)
{
    ConfigParser config;
    if (!config.parse(filename)) {
        std::cerr << "Error: Could not parse file " << filename << std::endl;
        return false;
    }

    auto parseDoubleRange = [](const std::string& str) -> DoubleRange {
        size_t pos = str.find('(');
        size_t posEnd = str.find(')');
        if (pos != std::string::npos && posEnd != std::string::npos) {
            std::string range = str.substr(pos + 1, posEnd - pos - 1);
            size_t commaPos = range.find(',');
            double lower = std::stod(range.substr(0, commaPos));
            double upper = std::stod(range.substr(commaPos + 1));
            return DoubleRange(lower, upper);
        }
        return DoubleRange(0.0, 0.0);
    };

    auto parseIntRange = [](const std::string& str) -> IntRange {
        size_t pos = str.find('(');
        size_t posEnd = str.find(')');
        if (pos != std::string::npos && posEnd != std::string::npos) {
            std::string range = str.substr(pos + 1, posEnd - pos - 1);
            size_t commaPos = range.find(',');
            int lower = std::stoi(range.substr(0, commaPos));
            int upper = std::stoi(range.substr(commaPos + 1));
            return IntRange(lower, upper);
        }
        return IntRange(0, 0);
    };
    

    // Database
    sDatabaseName = config.get("Database", "database_name");
    bDecoySearch = config.get("Database", "decoy_search") == "1";
    iDecoyKind = std::stoi(config.get("Database", "decoy_kind"));

    // Tolerance
    precursorTolerance = parseDoubleRange(config.get("Tolerance", "precursor_tolerance"));
    fragmentTolerance = parseDoubleRange(config.get("Tolerance", "fragment_tolerance"));
    iIsotopeError = std::stoi(config.get("Tolerance", "isotope_error"));
    bDeoxidation = config.get("Tolerance", "deoxidation") == "1";
    
    // Charge
    fragmentChargeRange = parseIntRange(config.get("Charge", "fragment_charge_range"));
    precursorChargeRange = parseIntRange(config.get("Charge", "precursor_charge_range"));
    bPositiveCharge = config.get("Charge", "positive_charge") == "1";
    sAdductIon = config.get("Charge", "adduct_ion");

    // Thread
    iThreadNum = std::stoi(config.get("Thread", "thread_num"));

    // Enzyme
    sEnzymeDatabase = config.get("Enzyme", "enzyme_database");
    sEnzymeName = config.get("Enzyme", "enzyme_name");
    iAllowedMissed = std::stoi(config.get("Enzyme", "allowed_missed_cleavage"));

    // Modification
    sVarModsDatabase = config.get("Modification", "variable_mods_database");
    sStaticModsDatabase = config.get("Modification", "static_mods_database");
    sOligo3Mod = config.get("Modification", "oligo_3_mod");
    sOligo5Mod = config.get("Modification", "oligo_5_mod");
    sNucleic3Mod = config.get("Modification", "nucleic_3_mod");
    sNucleic5Mod = config.get("Modification", "nucleic_5_mod");
    iMaxVarMods = std::stoi(config.get("Modification", "max_var_mods"));

    // Fragment
    bUseAIons = config.get("Fragment", "use_a_ions") == "1";  
    bUseBIons = config.get("Fragment", "use_b_ions") == "1";
    bUseCIons = config.get("Fragment", "use_c_ions") == "1";
    bUseDIons = config.get("Fragment", "use_d_ions") == "1";
    bUseWIons = config.get("Fragment", "use_w_ions") == "1";
    bUseXIons = config.get("Fragment", "use_x_ions") == "1";
    bUseYIons = config.get("Fragment", "use_y_ions") == "1";
    bUseZIons = config.get("Fragment", "use_z_ions") == "1";
    bUseABIons = config.get("Fragment", "use_a_b_ions") == "1";
    bUseNLIons = config.get("Fragment", "use_nl_ions") == "1";
    fragmentMassRange = parseDoubleRange(config.get("Fragment", "fragment_mass_range"));

    // Oligonucleotide 
    digestMassRange = parseDoubleRange(config.get("Oligonucleotide", "digest_mass_range"));
    oligoLengthRange = parseIntRange(config.get("Oligonucleotide", "oligo_length_range"));

    // Spectrum
    scanRange = parseIntRange(config.get("Spectrum", "scan_range"));
    dMinimumIntensity = std::stod(config.get("Spectrum", "minimum_intensity"));
    bRemovePrecursorPeak = config.get("Spectrum", "remove_precursor_peak") == "1";
    dRemovePrecursorTolerance = std::stod(config.get("Spectrum", "remove_precursor_tolerance"));
    clearMZRange = parseDoubleRange(config.get("Spectrum", "clear_mz_range"));
    dPercentageBasePeak = std::stod(config.get("Spectrum", "percentage_base_peak"));
    iMSLevel = std::stoi(config.get("Spectrum", "ms_level"));
    iSpectrumBatchSize = std::stoi(config.get("Spectrum", "spectrum_batch_size"));

    // Search
    dBinWidth = std::stod(config.get("Search", "bin_width"));
    iXcorrProcessingOffset = std::stoi(config.get("Search", "xcorr_offset"));

    return true;
}

Options& Options::operator=(const Options& other)
{
    if (this != &other) {
        sDatabaseName = other.sDatabaseName;
        bDecoySearch = other.bDecoySearch;
        iDecoyKind = other.iDecoyKind;

        precursorTolerance = other.precursorTolerance;
        fragmentTolerance = other.fragmentTolerance;
        iIsotopeError = other.iIsotopeError;
        bDeoxidation = other.bDeoxidation;

        fragmentChargeRange = other.fragmentChargeRange;
        precursorChargeRange = other.precursorChargeRange;
        bPositiveCharge = other.bPositiveCharge;
        sAdductIon = other.sAdductIon;

        iThreadNum = other.iThreadNum;

        sEnzymeDatabase = other.sEnzymeDatabase;
        sEnzymeName = other.sEnzymeName;
        iAllowedMissed = other.iAllowedMissed;

        sVarModsDatabase = other.sVarModsDatabase;
        sStaticModsDatabase = other.sStaticModsDatabase;
        sOligo3Mod = other.sOligo3Mod;
        sOligo5Mod = other.sOligo5Mod;
        sNucleic3Mod = other.sNucleic3Mod;
        sNucleic5Mod = other.sNucleic5Mod;
        iMaxVarMods = other.iMaxVarMods;

        bUseAIons = other.bUseAIons;
        bUseBIons = other.bUseBIons;
        bUseCIons = other.bUseCIons;
        bUseDIons = other.bUseDIons;
        bUseWIons = other.bUseWIons;
        bUseXIons = other.bUseXIons;
        bUseYIons = other.bUseYIons;
        bUseZIons = other.bUseZIons;
        bUseABIons = other.bUseABIons;
        bUseNLIons = other.bUseNLIons;

        digestMassRange = other.digestMassRange;
        fragmentMassRange = other.fragmentMassRange;
        oligoLengthRange = other.oligoLengthRange;

        scanRange = other.scanRange;
        dMinimumIntensity = other.dMinimumIntensity;
        bRemovePrecursorPeak = other.bRemovePrecursorPeak;
        dRemovePrecursorTolerance = other.dRemovePrecursorTolerance;
        clearMZRange = other.clearMZRange;
        dPercentageBasePeak = other.dPercentageBasePeak;
        iMSLevel = other.iMSLevel;
        iSpectrumBatchSize = other.iSpectrumBatchSize;
    }
    return *this;
}

bool Options::printDefaultOptions(const std::string& sParamFilename)
{
    std::ofstream paramFile(sParamFilename);
    if (!paramFile) {
        std::cerr << "Error: Cannot write to file " << sParamFilename << std::endl;
        return false;
    }

    paramFile << "[Version]\n";
    paramFile << "version = 2024.12_rev_0\n\n";

    paramFile << "[Database]\n";
    paramFile << "database_name = ./resource/xxx.fasta\n";
    paramFile << "decoy_search = 1                       # 0=no (default), 1=decoy\n";
    paramFile << "decoy_kind = 1                         # 0=reverse, 1=shuffle\n\n";

    paramFile << "[Tolerance]\n";
    paramFile << "precursor_tolerance = (-16.0, 16.0)                    # upper bound of the precursor mass tolerance\n";
    paramFile << "fragment_tolerance = (-16.0, 16.0)                     # upper bound of the precursor mass tolerance\n";
    paramFile << "isotope_error = 4                                   # 0=off, 1=0/1 (C13 error), 2=0/1/2, 3=0/1/2/3, 4=-1/0/1/2/3, 5=-1/0/1/2/3\n";
    paramFile << "deoxidation = 0                                     # 0 = RNA, 1 = DNA\n\n";

    paramFile << "[Charge]\n";
    paramFile << "fragment_charge_range = (1, 3)                 # set maximum fragment charge state to analyze (allowed max 5)\n";
    paramFile << "precursor_charge_range = (1, 6)                # set maximum precursor charge state to analyze (allowed max 9)\n";
    paramFile << "positive_charge = 0                            # 0 = negative, 1 = positive\n";
    paramFile << "adduct_ion = H-							       # choice: H+, Na+, NH4+, K+, H-, only one can be selected.\n\n";

    paramFile << "[Thread]\n";
    paramFile << "thread_num = 2                         # 0 = poll CPU to set num threads; else specify num threads directly (max 128)\n\n";

    paramFile << "[Enzyme]\n";
    paramFile << "enzyme_database = ./resource/enzyme.ini\n";
    paramFile << "enzyme_name = enzyme_A                 # choose from list at end of this params file\n";
    paramFile << "allowed_missed_cleavage = 2             # maximum value is 5; for enzyme search\n\n";

    paramFile << "[Modification]\n";
    paramFile << "variable_mods_database = ./resource/varMod.ini                   \n";
    paramFile << "static_mods_database = ./resource/staticMod.ini\n";
    paramFile << "oligo_3_mod = OH,cP                                              # option = [\"OH\", \"cP\"], can choose multiple\n";
    paramFile << "oligo_5_mod = OH,P                                              # option = [\"OH\", \"P\"]\n";
    paramFile << "nucleic_3_mod = OH,cP                                                 # option = [\"OH\", \"cP\"]\n";
    paramFile << "nucleic_5_mod = OH,P                                                 # option = [\"OH\", \"P\"]\n";
    paramFile << "max_var_mods = 5                                    # limit variable mods to subset of specified proteins if this file is specified & present\n\n";

    paramFile << "[Fragment]\n";
    paramFile << "use_a_ions = 1\n";
    paramFile << "use_b_ions = 1\n";
    paramFile << "use_c_ions = 1\n";
    paramFile << "use_d_ions = 1\n";
    paramFile << "use_w_ions = 1\n";
    paramFile << "use_x_ions = 1\n";
    paramFile << "use_y_ions = 1\n";
    paramFile << "use_z_ions = 1\n";
    paramFile << "use_a_B_ions = 1\n";
    paramFile << "use_NL_ions = 1                        # 0=no, 1=yes to consider NH3/H2O neutral loss peaks\n\n";

    paramFile << "[Oligonucleotide]\n";
    paramFile << "digest_mass_range = (200.0, 5000.0)                # MH+ peptide mass range to analyze \n";
    paramFile << "fragment_mass_range = (200.0, 5000.0)              # fragment mass range to analyze\n";
    paramFile << "oligo_length_range = (3, 50)             # minimum and maximum peptide length to analyze (default min 1 to allowed max 51)\n\n";

    paramFile << "[Spectrum]\n";
    paramFile << "scan_range = (0, 0)                       # start and end scan range to search; either entry can be set independently\n";
    paramFile << "minimum_intensity = 0                 # minimum intensity value to read in\n";
    paramFile << "remove_precursor_peak = 0              # 0=no, 1=yes, 2=all charge reduced precursor peaks (for ETD), 3=phosphate neutral loss peaks\n";
    paramFile << "remove_precursor_tolerance = 1.5       # +- Da tolerance for precursor removal\n";
    paramFile << "clear_mz_range = (0.0, 0.0)               # clear out all peaks in the specified m/z range e.g. remove reporter ion region of TMT spectra\n";
    paramFile << "percentage_base_peak = 0.0             # specify a percentage (e.g. \"0.05\" for 5%) of the base peak intensity as a minimum intensity threshold\n";
    paramFile << "ms_level = 2                           # MS level to analyze, valid are levels 2 (default) or 3\n";
    paramFile << "spectrum_batch_size = 15000            # max. # of spectra to search at a time; 0 to search the entire scan range in one loop\n";

    paramFile.close();
    return true;
}

// Param

Param::Param()
{
    initMutex();
}

Param::~Param()
{
    destroyMutex();
}

void Param::initMutex()
{
    Threading::CreateMutex(&mutex_vInputFile);

    Threading::CreateMutex(&mutex_vNucleicAcid);

    Threading::CreateMutex(&mutex_vOligonucleotide);

    Threading::CreateMutex(&mutex_vVarModFormulas);
}

void Param::destroyMutex()
{
    Threading::DestroyMutex(mutex_vInputFile);

    Threading::DestroyMutex(mutex_vNucleicAcid);

    Threading::DestroyMutex(mutex_vOligonucleotide);

    Threading::DestroyMutex(mutex_vVarModFormulas);
}

bool Param::loadOptions(const std::string& sParamFilename)
{
    return options.loadParameters(sParamFilename);
}

bool Param::printDefaultOptions(const std::string& sParamFilename)
{
    return options.printDefaultOptions(sParamFilename);
}

// bool Param::ProcessCmdLine(int argc, char* argv[])
// {
//     vector<string> vArgv;
//     int iPrintParams = 0;
//     int iArgvIndex = 1;
//     string sTmpStr;
//     string sParamsFilename;

//     if (argc == 1)
//     {
//         string sErrorMsg = "\n Nebulion version %s\n\n" + g_sNebulionVersion + " Error - no input files specified so nothing to do.\n";
//         logerr(sErrorMsg.c_str());
//         exit(1);
//     }

//     sParamsFilename = "./nebulion.params";                    //default parameter file name, may change in SetOptions(...) function

//     //read all args
//     while (iArgvIndex < argc)
//     {
//         sTmpStr = string(argv[iArgvIndex]);
//         vArgv.push_back(sTmpStr);
//         ++iArgvIndex;
//     }

//     //reload parameter file
//     for (string& arg : vArgv)
//     {
//         if (arg[0] == '-' && arg[1] == 'P')
//             sParamsFilename = arg.substr(2);
//         if (arg[0] == '-' && arg[1] == 'p')
//             Param::printParams("./nebulion.params.new", 1);
//         if (arg[0] == '-' && arg[1] == 'q')
//             Param::printParams("./nebulion.params.new", 2);
//     }

//     if (!Param::loadParam(sParamsFilename))
//     {
//         string sErrorMsg = "\n Nebulion version " + g_sNebulionVersion + "\n\n" + " Error - load parameter file fail. \n";
//         logerr(sErrorMsg.c_str());
//         return false;
//     }

//     //set options and read inputfile
//     for (string& arg : vArgv)
//     {
//         if (arg[0] == '-')
//             Param::SetOptions(arg);
//         else
//         {
//             if (!IsValidInputFile(arg))
//             {
//                 string sErrorMsg = "\n Nebulion version " + g_sNebulionVersion + "\n\n" + " Error - input file \"" + arg + "\" not found.\n";
//                 logerr(sErrorMsg.c_str());
//                 return false;
//             }
//             InputFile inputFile;
//             inputFile.sFileName = arg;
//             inputFile.iFirstScan = Param::g_staticParams.options.scanRange.iStart;
//             inputFile.iLastScan = Param::g_staticParams.options.scanRange.iEnd;
//             Param::g_vInputFiles.push_back(inputFile);
//         }
//     }

//     return true;
// }

// bool Param::printUsage(string sCmd)
// {
//     string sLogMsg;

//     logout("\n");
//     sLogMsg = " Nebulion version \" " + g_sNebulionVersion + "\"\n";
//     logout(sLogMsg.c_str());
//     logout("\n");
//     sLogMsg = " Nebulion usage:  " + sCmd + "\"\n";
//     logout(sLogMsg.c_str());
//     logout("\n");
//     logout(" Supported input formats include mzXML, mzML, Thermo raw, mgf, and ms2 variants (cms2, bms2, ms2)\n");

//     logout("\n");
//     logout("       options:  -p         to print out a nebulion.params.new file\n");
//     logout("                 -q         to print out a nebulion.params.new file with more parameter entries\n");
//     logout("                 -P<params> to specify an alternate parameters file (default nebulion.params)\n");
//     logout("                 -D<dbase>  to specify a sequence database, overriding entry in parameters file\n");
//     logout("                 -F<num>    to specify the first/start scan to search, overriding entry in parameters file\n");
//     logout("                 -L<num>    to specify the last/end scan to search, overriding entry in parameters file\n");
//     logout("                            (-L option is required if -F option is used)\n");
//     logout("\n");
//     sLogMsg = "       example:  " + sCmd + "file1.mzXML file2.mzXML\n";
//     logout(sLogMsg.c_str());
//     sLogMsg = "            or   " + sCmd + "-F1000 -L1500 file1.mzXML    <- to search scans 1000 through 1500\n";
//     logout(sLogMsg.c_str());
//     sLogMsg = "            or   " + sCmd + "-PParams.txt *.mzXML         <- use parameters in the file 'Params.txt'\n";
//     logout(sLogMsg.c_str());
//     logout("\n");

//     exit(1);
// }

// bool Param::SetOptions(string arg)
// {
//     switch (arg[1])
//     {
//     case 'D':                        // Alternate sequence database.
//         if (arg.size() <= 2)
//             logerr("Missing text for parameter option -D<database>.  Ignored.\n");
//         else
//             Param::g_staticParams.databaseInfo.sFilename = arg.substr(2);
//         break;
//     case 'F':   // first scan
//         if (arg.size() <= 2)
//             logerr("Missing text for parameter option -F<num>.  Ignored.\n");
//         else
//             Param::g_staticParams.options.scanRange.iStart = stoi(arg.substr(2));
//         break;
//     case 'L':  // last scan
//         if (arg.size() <= 2)
//             logerr("Missing text for parameter option -L<num>.  Ignored.\n");
//         else
//             Param::g_staticParams.options.scanRange.iEnd = stoi(arg.substr(2));
//         break;
//     case 'B':    // batch size
//         if (arg.size() <= 2)
//             logerr("Missing text for parameter option -B<num>.  Ignored.\n");
//         else
//             Param::g_staticParams.options.iSpectrumBatchSize = stoi(arg.substr(2));
//         break;
//     default:
//         break;
//     }
//     return true;
// }













