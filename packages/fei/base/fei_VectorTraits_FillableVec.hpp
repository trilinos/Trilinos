/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_VectorTraits_FillableVec_hpp_
#define _fei_VectorTraits_FillableVec_hpp_

#include <fei_VectorTraits.hpp>
#include <fei_FillableVec.hpp>

namespace fei {
  template<>
  struct VectorTraits<FillableVec> {
    static const char* typeName()
      { return("FillableVec"); }
    
    static int setValues(FillableVec* vec, int firstLocalOffset,
                         double scalar, bool isSolnVector=false)
      {
        vec->setValues(scalar);
        return(0);
      }

    static int putValuesIn(FillableVec* vec,
                     int firstLocalOffset,
                     int numValues, const int* indices, const double* values,
                     bool sum_into,
                     bool isSolnVector=false,
                     int vectorIndex=0)
      {
        if (sum_into) {
          for(int i=0; i<numValues; ++i) {
            if (indices[i] < 0) continue;
            vec->addEntry(indices[i], values[i]);
          }
        }
        else {
          for(int i=0; i<numValues; ++i) {
            if (indices[i] < 0) continue;
            vec->putEntry(indices[i], values[i]);
          }
        }

        return( 0 );
      }

    static int copyOut(FillableVec* vec,
                       int firstLocalOffset,
                       int numValues, const int* indices, double* values,
                       bool isSolnVector=false,
                       int vectorIndex=0)
      {
        for(int i=0; i<numValues; ++i) {
          try {
            values[i] = vec->getEntry(indices[i]);
          }
          catch(...) {}
        }
        return(0);
      }

    static int update(FillableVec* vec,
                      double a,
                      FillableVec* x,
                      double b)
    { return(-1); }

  };
}//namespace fei

#endif // _fei_VectorTraits_FillableVec_hpp_
