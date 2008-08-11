#ifndef _BCNodeSet_h_
#define _BCNodeSet_h_

/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <fei_macros.hpp>

#include <fei_defs.h>

class BCNodeSet {
 public:
   BCNodeSet();
   ~BCNodeSet();

   int numNodes_;
   GlobalID* nodeIDs_;
   int fieldID_;
   int* offsetsIntoField_;
   double* prescribed_values_;

 private:
   void deleteMemory();
};

#endif

