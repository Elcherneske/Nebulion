// Copyight 2024 Yucheng Liao
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


#include "Common.h"
#include "Param.h"
#include "NebulionInterface.h"

#include <iostream>
#include <string>

using namespace NebulionInterface;

bool ProcessCmdLine(int argc, char* argv[]);

int main(int argc, char* argv[])
{
    bool bSearchSucceeded = false;

    ISearchManager* pSearchMgr = GetSearchManager();

    Param::getInstance().loadOptions("D:/code/VSC++/Nebulion/resource/params.config");

    ProcessCmdLine(argc, argv);  

    bSearchSucceeded = pSearchMgr->DoSearch();

    NebulionInterface::ReleaseSearchManager();

    if (!bSearchSucceeded)
    {
        exit(1);
    }

    return (0);
} // main

bool ProcessCmdLine(int argc, char* argv[])
{
    if (argc == 1)
    {
        std::cerr << " Error - no input files specified so nothing to do.\n";
        exit(1);
    }

    //set options and read inputfile
    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];
        InputFile inputFile;
        inputFile.sFileName = arg;
        Param::getInstance().vInputFiles.push_back(inputFile);
        
    }

    return true;
}

