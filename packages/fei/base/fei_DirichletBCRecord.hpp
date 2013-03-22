#ifndef _fei_DirichletBCRecord_hpp_
#define _fei_DirichletBCRecord_hpp_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

namespace fei {

struct DirichletBCRecord {
  int IDType;
  int ID;
  int fieldID;
  int whichComponent;
  double prescribedValue;

  bool operator!=(const DirichletBCRecord& rhs) const
  {
    return IDType != rhs.IDType || ID != rhs.ID || 
           fieldID != rhs.fieldID || whichComponent != rhs.whichComponent;
  }
};

class less_DirichletBCRecord {
 public:
  less_DirichletBCRecord(){}
  ~less_DirichletBCRecord(){}

  bool operator()(const DirichletBCRecord& lhs,
                  const DirichletBCRecord& rhs)
  {
    if (lhs.IDType < rhs.IDType) return true;
    if (lhs.IDType > rhs.IDType) return false;

    if (lhs.ID < rhs.ID) return true;
    if (lhs.ID > rhs.ID) return false;

    if (lhs.fieldID < rhs.fieldID) return true;
    if (lhs.fieldID > rhs.fieldID) return false;

    if (lhs.whichComponent < rhs.whichComponent) return true;
    if (lhs.whichComponent > rhs.whichComponent) return false;

    return false;
  }
};

}//namespace fei

#endif

