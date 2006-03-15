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
