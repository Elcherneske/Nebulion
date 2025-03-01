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


#ifndef _NEBULIONINTERFACES_H_
#define _NEBULIONINTERFACES_H_

#include "Common.h"
#include "DataStructure.h"
#include "ThreadPool.h"

namespace NebulionInterface
{
   class ISearchManager
   {
   public:
      virtual ~ISearchManager() {}
      virtual bool DoSearch() = 0;
      virtual bool IsSearchError() = 0;
      virtual void GetStatusMessage(string& strStatusMsg) = 0;
      virtual void CancelSearch() = 0;
      virtual bool IsCancelSearch() = 0;
      virtual void ResetSearchStatus() = 0;
   };

   ISearchManager* GetSearchManager();
   void ReleaseSearchManager();
   static ThreadPool* _tp;
}

#endif
