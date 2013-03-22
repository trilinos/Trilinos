/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_ReverseMapper_hpp_
#define _fei_ReverseMapper_hpp_

#include <fei_macros.hpp>

#include <map>

#include <fei_EqnRecord.hpp>

namespace fei {
  class VectorSpace;

/** Allows mapping from equation-numbers to IDs, fields, etc. */
class ReverseMapper {
 public:
  /** constructor */
  ReverseMapper(const VectorSpace& vspace);

  /** destructor */
  virtual ~ReverseMapper();

  EqnRecord getEqnRecord(int global_eqn, int option=0) const;

 private:
  std::map<int,EqnRecord> eqnmap_;

  ReverseMapper(const ReverseMapper& src);
  ReverseMapper& operator=(const ReverseMapper& src);
};//class ReverseMapper
}//namespace fei
#endif

