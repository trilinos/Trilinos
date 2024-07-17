// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
  Exception(const std::string FileName, const int LineNumber,
            const std::string Line1, const std::string Line2 = "",
            const std::string Line3 = "", const std::string Line4 = "",
            const std::string Line5 = "", const std::string Line6 = "") :
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
  std::string FileName_;
  int LineNumber_;
  std::string Line1_;
  std::string Line2_;
  std::string Line3_;
  std::string Line4_;
  std::string Line5_;
  std::string Line6_;

}; // class Exception

} // namespace Galeri

#endif
