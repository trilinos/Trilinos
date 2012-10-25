/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_VectorTraits_FEData_hpp_
#define _fei_VectorTraits_FEData_hpp_

//This file defines vector traits for FiniteElementData vectors
//(well, "vector-views" to be more precise).
//

#include <fei_VectorTraits.hpp>
#include <fei_FiniteElementData.hpp>

namespace fei {

  /** specialization for FiniteElementData */
  template<>
  struct VectorTraits<FiniteElementData>  {

    /** name of VectorTraits type */
    static const char* typeName()
      { return("FiniteElementData"); }

    /** set all vector values to specified scalar */
    static int setValues(FiniteElementData* vec, int firstLocalOffset,
			 double scalar, bool isSolnVector=false)
      {
	return(-1);
      }

    /** sum-into operation for vector data */
    static int putValuesIn(FiniteElementData* vec,
		     int firstLocalOffset,
		     int numValues, const int* indices, const double* values,
                     bool sum_into,
		     bool isSolnVector=false,
		     int vectorIndex=0)
      {
	return(-1);
      }

    /** copy out vector data */
    static int copyOut(FiniteElementData* vec,
		       int firstLocalOffset,
		       int numValues, const int* indices, double* values,
		       bool isSolnVector=false,
		       int vectorIndex=0)
      {
	return(-1);
      }

    /** vec = b*vec + a*x */
    static int update(FiniteElementData* vec,
		      double a,
		      const FiniteElementData* x,
		      double b)
    { return(-1); }

  };//struct VectorTraits
}//namespace fei

#endif // _fei_VectorTraits_FEData_hpp_
