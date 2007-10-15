/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_fwd_hpp_
#define _snl_fei_fwd_hpp_

#include <fei_fwd.hpp>

namespace snl_fei {
  template<class T> class CommUtils; 
  template<class RecordType,class RecordType_COMPARE=lessthan<int> > class Constraint;
  class RecordCollection;

  class BlockDescriptor;
  class PointBlockMap;

  class Broker;
  class Factory;
} //namespace snl_fei

#endif // _snl_fei_fwd_hpp_
