/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_VectorTraits_CSVec_hpp_
#define _fei_VectorTraits_CSVec_hpp_

#include <fei_VectorTraits.hpp>
#include <fei_CSVec.hpp>

namespace fei {
  template<>
  struct VectorTraits<CSVec> {
    static const char* typeName()
      { return("CSVec"); }
    
    static int setValues(CSVec* vec, int firstLocalOffset,
                         double scalar, bool isSolnVector=false)
      {
        set_values(*vec, scalar);
        return(0);
      }

    static int putValuesIn(CSVec* vec,
                     int firstLocalOffset,
                     int numValues, const int* indices, const double* values,
                     bool sum_into,
                     bool isSolnVector=false,
                     int vectorIndex=0)
      {
        if (sum_into) {
          for(int i=0; i<numValues; ++i) {
            if (indices[i] < 0) continue;
            add_entry(*vec, indices[i], values[i]);
          }
        }
        else {
          for(int i=0; i<numValues; ++i) {
            if (indices[i] < 0) continue;
            put_entry(*vec, indices[i], values[i]);
          }
        }

        return( 0 );
      }

    static int copyOut(CSVec* vec,
                       int firstLocalOffset,
                       int numValues, const int* indices, double* values,
                       bool isSolnVector=false,
                       int vectorIndex=0)
      {
        for(int i=0; i<numValues; ++i) {
          try {
            values[i] = get_entry(*vec, indices[i]);
          }
          catch(...) {}
        }
        return(0);
      }

    static int update(CSVec* vec,
                      double a,
                      const CSVec* x,
                      double b)
    { return(-1); }

  };
}//namespace fei

#endif // _fei_VectorTraits_CSVec_hpp_

