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


#ifndef _SPECTRUMPREPROCESS_H_
#define _SPECTRUMPREPROCESS_H_

#include "Common.h"
#include "DataStructure.h"
#include "Param.h"

class SpectrumPreprocessUtils
{
public:
    SpectrumPreprocessUtils();
   ~SpectrumPreprocessUtils();

   static MSReader& getMSReader();

   static bool LoadAllSpectrum(vector<InputFile>& inputFiles, vector<std::pair<InputFile, vector<CometSpectrum>>>& vInputFileCometSpectrumPair);

   static bool LoadSpectrum(InputFile& inputFile, vector<CometSpectrum>& vCometSpectrumList);

   static bool PreprocessSpectrum(Spectrum& spec, vector<CometSpectrum>& cometSpectrumList);

private:

    static MSReader msreader;

    static bool CalculateMassTolerance(CometSpectrum& cometSpectrum);

};

#endif // _COMETSPECTRUMPREPROCESS_H_