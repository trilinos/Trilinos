/*
//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

#ifndef EPETRAEXT_EXCEPTION
#define EPETRAEXT_EXCEPTION

#include "EpetraExt_ConfigDefs.h"
#include <string>

namespace EpetraExt
{

class Exception 
{
public:
  Exception(const string FileName, const int LineNumber,
            const string Line1, const string Line2 = "",
            const string Line3 = "", const string Line4 = "",
            const string Line5 = "", const string Line6 = "") :
    FileName_(FileName),
    LineNumber_(LineNumber),
    Line1_(Line1),
    Line2_(Line2),
    Line3_(Line3),
    Line4_(Line4),
    Line5_(Line5),
    Line6_(Line6)
  {}

  void Print()
  {
    cerr << "EpetraExt: Exception: " << Line1_ << endl;
    if (Line2_ != "")
      cerr << "EpetraExt: Exception: " << Line2_ << endl;
    if (Line3_ != "")
      cerr << "EpetraExt: Exception: " << Line3_ << endl;
    if (Line4_ != "")
      cerr << "EpetraExt: Exception: " << Line4_ << endl;
    if (Line5_ != "")
      cerr << "EpetraExt: Exception: " << Line5_ << endl;
    if (Line6_ != "")
      cerr << "EpetraExt: Exception: " << Line6_ << endl;

    cerr << "EpetraExt: File: " << FileName_ << ", line: " << LineNumber_ << endl;
  }

private:
  string FileName_;
  int LineNumber_;
  string Line1_;
  string Line2_;
  string Line3_;
  string Line4_;
  string Line5_;
  string Line6_;

}; // class Exception

} // namespace EpetraExt

#endif
