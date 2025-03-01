// Copyright 2023 Jimmy Eng
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


#ifndef _STATUS_H_
#define _STATUS_H_

#include "Common.h"
#include "Threading.h"

enum Result
{
    Result_Succeeded = 0,
    Result_Failed,
    Result_Cancelled
};

class Status
{
public:
    Status()
   {
       _cometResult = Result_Succeeded;
      _strStatusMessage = "";

      Threading::CreateMutex(&_statusCheckMutex);
   }

   ~Status()
   {
      Threading::DestroyMutex(_statusCheckMutex);
   }

   Result GetStatus()
   {
      Result result;
      Threading::LockMutex(_statusCheckMutex);
      result = _cometResult;
      Threading::UnlockMutex(_statusCheckMutex);

      return result;
   }

   Result GetStatus(std::string& strStatusMsg)
   {
       Result result;
      Threading::LockMutex(_statusCheckMutex);
      result = _cometResult;
      strStatusMsg = _strStatusMessage;
      Threading::UnlockMutex(_statusCheckMutex);

      return result;
   }

   void SetStatus(Result result)
   {
      Threading::LockMutex(_statusCheckMutex);
      _cometResult = result;
      Threading::UnlockMutex(_statusCheckMutex);
   }

   void SetStatus(Result result, const std::string &strStatusMsg)
   {
      Threading::LockMutex(_statusCheckMutex);
      _cometResult = result;
      _strStatusMessage = strStatusMsg;
      Threading::UnlockMutex(_statusCheckMutex);
   }

   void GetStatusMsg(std::string& strStatusMsg)
   {
      Threading::LockMutex(_statusCheckMutex);
      strStatusMsg = _strStatusMessage;
      Threading::UnlockMutex(_statusCheckMutex);
   }

   void SetStatusMsg(const std::string &strStatusMsg)
   {
      Threading::LockMutex(_statusCheckMutex);
      _strStatusMessage = strStatusMsg;
      Threading::UnlockMutex(_statusCheckMutex);
   }

   void ResetStatus()
   {
      Threading::LockMutex(_statusCheckMutex);
      _cometResult = Result_Succeeded;
      _strStatusMessage = "";
      Threading::UnlockMutex(_statusCheckMutex);
   }

   bool IsError()
   {
      bool bError;
      Threading::LockMutex(_statusCheckMutex);
      bError = _cometResult == Result_Failed;
      Threading::UnlockMutex(_statusCheckMutex);

      return bError;
   }

   bool IsCancel()
   {
      bool bCancelled;
      Threading::LockMutex(_statusCheckMutex);
      bCancelled = _cometResult == Result_Cancelled;
      Threading::UnlockMutex(_statusCheckMutex);

      return bCancelled;
   }

private:
    Result _cometResult;
   std::string _strStatusMessage;
   Mutex _statusCheckMutex;
};

#endif // _COMETSTATUS_H_
