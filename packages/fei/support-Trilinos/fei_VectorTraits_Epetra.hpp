/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_VectorTraits_Epetra_h_
#define _fei_VectorTraits_Epetra_h_

//
//IMPORTANT NOTE: Make sure that wherever this file is included from, it
//appears BEFORE any include of fei_base.hpp or fei_Vector.hpp !!!
//
#include <fei_VectorTraits.hpp>
#include <fei_Include_Trilinos.hpp>

namespace fei {
  /** Declare an Epetra_MultiVector specialization of the
      fei::VectorTraits struct.

      This allows Epetra_MultiVector to be used as the template parameter of the
      fei::Vector class.
  */
  template<>
  struct VectorTraits<Epetra_MultiVector> {
    static const char* typeName()
      { return("Epetra_MultiVector"); }

    static int setValues(Epetra_MultiVector* vec, int firstLocalOffset,
			 double scalar, bool isSolnVector=false)
      {
	return( vec->PutScalar(scalar) );
      }

    static int putValuesIn(Epetra_MultiVector* vec,
		           int firstLocalOffset,
		           int numValues,
                           const int* indices,
                           const double* values,
                           bool sum_into,
                           bool isSolnVector=false,
                           int vectorIndex=0)
      {
        if (sum_into) {
          for(int i=0; i<numValues; ++i) {
            if (vec->SumIntoGlobalValue(indices[i],vectorIndex,values[i])!=0){
              return(-1);
            }
          }
	}
        else {
	  for(int i=0; i<numValues; ++i) {
	    if (vec->ReplaceGlobalValue(indices[i],vectorIndex,values[i])!=0){
	      return(-1);
	    }
	  }
        }
	return(0);
      }

    static int copyOut(Epetra_MultiVector* vec,
		       int firstLocalOffset,
		       int numValues, const int* indices, double* values,
		       bool isSolnVector=false,
		       int vectorIndex=0)
      {
	const Epetra_BlockMap& map = vec->Map();
	double* coefPtr = (*vec)[vectorIndex];
	for(int i=0; i<numValues; ++i) {
	  values[i] = coefPtr[map.LID(indices[i])];
	}

	return(0);
      }

    static double* getLocalCoefsPtr(Epetra_MultiVector* vec,
                                    bool isSolnVector=false,
                                    int vectorIndex=0)
      {
        return(vec->Pointers()[vectorIndex]);
      }

    static int update(Epetra_MultiVector* vec,
		      double a,
		      Epetra_MultiVector* x,
		      double b)
    {
      return( vec->Update(a, *x, b) );
    }

  };//struct VectorTraits<Epetra_MultiVector>
}//namespace fei

#endif // _fei_VectorTraits_Epetra_hpp_
