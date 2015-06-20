#ifndef geometry_transformH
#define geometry_transformH

#include "../rtcompiler/RTC_FunctionRTC.hh"
namespace PAMGEN_NEVADA {

class Geometry_Transform{
public:
  Geometry_Transform(const std::string & funcBody,
			   std::stringstream & error_stream);
  ~Geometry_Transform();
  void Display_Class(std::ostream&, const std::string &indent);
  void Operate(double* coords, long long num_nodes,long long dim);

  
private:
  std::string                     _funcBody;
  std::string                     application_direction;
  PG_RuntimeCompiler::Function  _function;

};
}
#endif
