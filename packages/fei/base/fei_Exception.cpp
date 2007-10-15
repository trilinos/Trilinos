/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_macros.hpp"
#include "fei_Exception.hpp"
#include "fei_iostream.hpp"
#include <cstdlib>

bool fei_exception_abort = false;

namespace fei {
Exception::Exception(const char* msg) throw()
  : std::exception(),
    smsg_(msg)
{
  if (fei_exception_abort) {
    FEI_CERR << "fei::Exception with msg: "<<smsg_<<
	 "User has requested fei::Exception to abort instead of throw."<<FEI_ENDL;
    std::abort();
  }
}

Exception::Exception(std::string msg) throw()
  : std::exception(),
    smsg_(msg)
{
  if (fei_exception_abort) {
    FEI_CERR << "fei::Exception with msg: "<<smsg_<<
	 "User has requested fei::Exception to abort instead of throw."<<FEI_ENDL;
    std::abort();
  }
}

Exception::~Exception() throw()
{
}

const char* Exception::what() const throw()
{
  return( smsg_.c_str() );
}

}//namespace fei

