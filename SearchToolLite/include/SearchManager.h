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


#ifndef _SEARCHMANAGER_H_
#define _SEARCHMANAGER_H_

#include "Common.h"
#include "DataStructure.h"
#include "Param.h"
#include "NebulionInterface.h"

using namespace NebulionInterface;

class SearchManager : public ISearchManager
{
public:
   SearchManager();
   ~SearchManager();

   // Methods inherited from ISearchManager
   virtual bool DoSearch();
   virtual void SetOutputFileBaseName(string pszBaseName);
   virtual bool IsSearchError();
   virtual void GetStatusMessage(string &strStatusMsg);
   virtual void CancelSearch();
   virtual bool IsCancelSearch();
   virtual void ResetSearchStatus();

};

#endif