// @HEADER
// ***********************************************************************
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
// @HEADER

#include <Teuchos_ParameterList.hpp>
#include <sstream>
#include <fstream>
#include <ostream>

namespace Z2{

void getOutputStreamFromParameterList(
  Teuchos::ParameterList &pl, std::string key, std::ostream *&os,
  std::ostream &defaultValue);

}  //namespace Z2
