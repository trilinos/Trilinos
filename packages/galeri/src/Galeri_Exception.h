// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

#ifndef GALERI_EXCEPTION
#define GALERI_EXCEPTION

#include "Galeri_ConfigDefs.h"
#include <string>

namespace Galeri
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
    cerr << "Galeri: Exception: " << Line1_ << endl;
    if (Line2_ != "")
      cerr << "Galeri: Exception: " << Line2_ << endl;
    if (Line3_ != "")
      cerr << "Galeri: Exception: " << Line3_ << endl;
    if (Line4_ != "")
      cerr << "Galeri: Exception: " << Line4_ << endl;
    if (Line5_ != "")
      cerr << "Galeri: Exception: " << Line5_ << endl;
    if (Line6_ != "")
      cerr << "Galeri: Exception: " << Line6_ << endl;

    cerr << "Galeri: File: " << FileName_ << ", line: " << LineNumber_ << endl;
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

} // namespace Galeri

#endif
