// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
