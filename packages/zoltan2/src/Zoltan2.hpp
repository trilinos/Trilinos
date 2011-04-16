
#include <string>
#include <sstream>
#include <iostream>

// Zoltan2 exceptions.  At the current time the sub classes add
//  nothing more, but they might in the future.

namespace Zoltan2{
  class Zoltan2_error{
    protected:
      std::string info;
      std::string fileName;
      int lineNumber;

    public:
      Zoltan2_error():info(),fileName(),lineNumber(0){}
      Zoltan2_error(const std::string i, const std::string f, 
        int n): info(i),fileName(f),lineNumber(n){}

      std::string what() const{
        std::ostringstream os;
        os << lineNumber;
        std::string s;
        if (fileName.size() > 0){
          s = s + fileName+" ("+os.str()+ ") ";
        }
        if (info.size() > 0){
          s = s + info;
        }
        return s;
      }
  };
  class memory_error: public Zoltan2_error {
    public:
      memory_error() {}
      memory_error(const std::string &i, const std::string &f, 
        int n)
      {
         info = i;
         fileName = f;
         lineNumber= n; 
      }
  };
  class usage_error: public Zoltan2_error {
    public:
      usage_error(){}
      usage_error(const std::string &i, const std::string &f, 
        int n)
      {
         info = i;
         fileName = f;
         lineNumber= n; 
      }
  };
  class data_type_error: public Zoltan2_error {
    public:
      data_type_error(){}
      data_type_error(const std::string &i, const std::string &f, 
        int n)
      {
         info = i;
         fileName = f;
         lineNumber= n; 
      }
  };
  class possible_bug: public Zoltan2_error {
    public:
      possible_bug(){}
      possible_bug(const std::string &i, const std::string &f, 
        int n)
      {
         info = i;
         fileName = f;
         lineNumber= n; 
      }
  };
  class unknown_error_type: public Zoltan2_error {
    public:
      unknown_error_type(){}
      unknown_error_type(const std::string &i, const std::string &f, 
        int n)
      {
         info = i;
         fileName = f;
         lineNumber= n; 
      }
  };
}

