/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_VectorTraits_LinProbMgr_hpp_
#define _fei_VectorTraits_LinProbMgr_hpp_


#include <fei_VectorTraits.hpp>
#include <fei_LinearProblemManager.hpp>

namespace fei {

/** This struct specialization defines vector traits for LinearProblemManager
  vector representations.
*/
  template<>
  struct VectorTraits<fei::LinearProblemManager>  {

   /** Return a string type-name for the vector. */
    static const char* typeName()
      { return("fei::LinearProblemManager"); }

   /** Set a specified scalar value throughout the vector.
       */
    static int setValues(fei::LinearProblemManager* vec, int firstLocalOffset,
                         double scalar, bool isSolnVector=false)
      {
        vec->setVectorValues(scalar, isSolnVector);
        return(0);
      }

   /** Sum values into the vector, adding to any
          that may already exist at the specified indices.
      */
    static int putValuesIn(fei::LinearProblemManager* vec,
		     int firstLocalOffset,
		     int numValues, const int* indices, const double* values,
                     bool sum_into,
		     bool isSolnVector=false,
		     int vectorIndex=0)
      {
        int err = vec->insertVectorValues(numValues, indices, values,
                                          sum_into, isSolnVector, vectorIndex);
	return(err);
      }

   /** Copy values from the specified indices out into the user-allocated
          array 'values'.
      */
    static int copyOut(fei::LinearProblemManager* vec,
		       int firstLocalOffset,
		       int numValues, const int* indices, double* values,
		       bool isSolnVector=false,
		       int vectorIndex=0)
      {
        int err = vec->copyOutVectorValues(numValues, indices, values,
                                           isSolnVector, vectorIndex);
        return(err);
      }

    /** Perform global assembly. */
    static int globalAssemble(fei::LinearProblemManager* vec)
      { return( vec->globalAssemble() ); }

    /** Update 'vec' = b*'vec' + a*x
       */
    static int update(fei::LinearProblemManager* vec,
		      double a,
		      const fei::LinearProblemManager* x,
		      double b)
    { return(-1); }

  };//struct VectorTraits
}//namespace fei

#endif // _fei_VectorTraits_LinProbMgr_hpp_
