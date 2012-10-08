/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_EqnRecord_hpp_
#define _fei_EqnRecord_hpp_

#include "fei_macros.hpp"

namespace fei {
/** Associate an equation with the IDType, ID and fieldID of the
  mesh-object it corresponds to. */
struct EqnRecord {
  /** Global equation index. */
  int global_eqn;

  /** IDType (usually corresponds to node, edge, element, etc.) */
  int IDType;

  /** Global ID of the node/edge/element/etc. */
  int ID;

  /** Identifier of the field that the global_eqn corresponds to. */
  int fieldID;

  /** Offset into the field. */
  int offset;
};//struct EqnRecord
}//namespace fei
#endif

