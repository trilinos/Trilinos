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
