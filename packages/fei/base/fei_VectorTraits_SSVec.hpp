/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_VectorTraits_SSVec_hpp_
#define _fei_VectorTraits_SSVec_hpp_

#include <fei_VectorTraits.hpp>

#include <fei_SSVec.hpp>

namespace fei {
  template<>
  struct VectorTraits<SSVec> {
    static const char* typeName()
      { return("SSVec"); }
    
    static int setValues(SSVec* vec, int firstLocalOffset,
			 double scalar, bool isSolnVector=false)
      {
	vec->coefs() = scalar;
	return(0);
      }

    static int putValuesIn(SSVec* vec,
		     int firstLocalOffset,
		     int numValues, const int* indices, const double* values,
                     bool sum_into,
		     bool isSolnVector=false,
		     int vectorIndex=0)
      {
	int err = 0;
        if (sum_into) {
          for(int i=0; i<numValues; ++i) {
            if (indices[i] < 0) continue;
            err = vec->addEntry(indices[i], values[i]);
            if (err) return(-1);
          }
        }
        else {
          for(int i=0; i<numValues; ++i) {
            if (indices[i] < 0) continue;
            err = vec->putEntry(indices[i], values[i]);
            if (err) return(-1);
          }
        }

	return(err);
      }

    static int copyOut(SSVec* vec,
		       int firstLocalOffset,
		       int numValues, const int* indices, double* values,
		       bool isSolnVector=false,
		       int vectorIndex=0)
      {
	for(int i=0; i<numValues; ++i) {
	  int index = vec->indices().find(indices[i]);
	  if (index < 0) continue;
	  values[i] = vec->coefs()[index];
	}
	return(0);
      }

    static int update(SSVec* vec,
		      double a,
		      SSVec* x,
		      double b)
    { return(-1); }

  };
}//namespace fei

#endif // _fei_VectorTraits_SSVec_hpp_
