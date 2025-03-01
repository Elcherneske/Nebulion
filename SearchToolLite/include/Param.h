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

#include <map>
#include <string>
#include <vector>

class ConfigParser {
private:
    std::map<std::string, std::map<std::string, std::string>> config;
public:
    bool parse(const std::string& filename);
    std::string get(const std::string& section, const std::string& option);
};

class TSVParser {
private:
    std::vector<std::map<std::string, std::string>> data;
    std::vector<std::string> headers;

    void parseHeaders(const std::string& headerLine);
    std::map<std::string, std::string> parseLine(const std::string& line);
public:
    bool readFile(const std::string& filename);

    std::vector<std::map<std::string, std::string>> getData() const;
    std::vector<std::string> getHeaders() const;
};

class Options
{
public:
   // Database params
   std::string sDatabaseName;
   bool bDecoySearch;
   int iDecoyKind;

   // Tolerance params
   DoubleRange precursorTolerance;
   DoubleRange fragmentTolerance; 
   int iIsotopeError;
   bool bDeoxidation;

   // Charge params
   IntRange fragmentChargeRange;
   IntRange precursorChargeRange;
   bool bPositiveCharge;
   std::string sAdductIon;

   // Thread params
   int iThreadNum;

   // Enzyme params
   std::string sEnzymeDatabase;
   std::string sEnzymeName;
   int iAllowedMissed;

   // Modification params
   std::string sVarModsDatabase;
   std::string sStaticModsDatabase;
   std::string sOligo3Mod;
   std::string sOligo5Mod;
   std::string sNucleic3Mod;
   std::string sNucleic5Mod;
   int iMaxVarMods;

   // Fragment params
   bool bUseAIons;
   bool bUseBIons;
   bool bUseCIons;
   bool bUseDIons;
   bool bUseWIons;
   bool bUseXIons;
   bool bUseYIons;
   bool bUseZIons;
   bool bUseABIons;
   bool bUseNLIons;
   bool bFreeRadical;

   // Oligonucleotide params
   DoubleRange digestMassRange; // the mass range of the digest oligonucleotide without var mods
   DoubleRange fragmentMassRange;
   IntRange oligoLengthRange;

   // Spectrum params
   IntRange scanRange;
   double dMinimumIntensity;
   bool bRemovePrecursorPeak;
   double dRemovePrecursorTolerance;
   DoubleRange clearMZRange;
   double dPercentageBasePeak;
   int iMSLevel;
   int iSpectrumBatchSize;

   // Search params
   double dBinWidth;
   int iXcorrProcessingOffset;
   Options();

   bool loadParameters(const std::string& filename);

   bool printDefaultOptions(const std::string& sParamFilename);

   Options& operator=(const Options& other);
};


class Param
{
public:
    Options                                      options;//user defined options
    std::string                                       sVersion;
    std::vector<InputFile>                            vInputFiles;
    std::vector<NucleicAcid>                          vNucleicAcids;
    std::vector<Oligonucleotide>                      vOligos;
    std::vector<VarModFormula>                          vVarModFormulas;
    Mutex                                        mutex_vInputFile;
    Mutex                                        mutex_vNucleicAcid;
    Mutex                                        mutex_vOligonucleotide;
    Mutex                                        mutex_vVarModFormulas;

    Param();   
    ~Param();

    bool loadOptions(const std::string& sParamFilename);
    bool printDefaultOptions(const std::string& sParamFilename);

    static Param& getInstance() {
        static Param param;
        return param;
    }

private:
    void initMutex();
    void destroyMutex();
};



#endif
