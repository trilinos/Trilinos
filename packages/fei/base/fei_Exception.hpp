/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_Exception_hpp_
#define _fei_Exception_hpp_

#include "fei_macros.hpp"
#include <string>
#include <exception>

namespace fei {

/** fei specialization of std::exception */
class Exception : public std::exception {
 public:
  /** constructor */
  Exception(const char* msg) throw();

  /** constructor */
  Exception(std::string msg) throw();

  /** destructor */
  virtual ~Exception() throw();

  /** return const char-ptr of exception message */
  const char* what() const throw();

 private:
  std::string smsg_;
};//class Exception

}//namespace fei

#endif // _fei_Exception_hpp_
