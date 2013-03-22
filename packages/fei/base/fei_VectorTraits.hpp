/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _fei_VectorTraits_hpp_
#define _fei_VectorTraits_hpp_

#include <fei_macros.hpp>

namespace fei {

  /** Define a struct of vector access traits. The fei vector implementation
      class fei::Vector is essentially a filter which passes data to
      library-specific vector objects (such as Trilinos/Epetra's
      Epetra_MultiVector).
      fei::Vector is a template, and the template parameter is the vector
      object. In order to use an arbitrary vector object with fei::Vector,
      it is necessary to define a specialization of this VectorTraits
      struct for the vector object.

      For an example specialization, see
      support-Trilinos/VectorTraits_Epetra.hpp.

      This "base" VectorTraits struct provides function stubs for default
      type "T", which will catch the use of any vector type for which
      specialized traits have not been defined.

      Note:
      Several functions have an argument called 'isSolnVector'. This allows
      for cases where the vector object is an aggregate container such as
      the old LinearSystemCore, or the newer LinearProblemManager, which
      embody a linear-system with 2 vectors (a solution vector and a right-
      hand-side vector).

      Note2:
      Most functions have an argument called 'firstLocalOffset'. This is
      the local processor's starting offset into the global index space.
      Other 'indices' arguments always contain global indices, so offsets
      into a processor's data can be obtained by, for example,
        'indices[i] - firstLocalOffset'.

      Note3:
      All functions are expected to only operate on the locally-stored
      portion of the vector. No inter-process communication should occur
      in any function except globalAssemble().
  */
  template<typename T>
    struct VectorTraits {

      /** Return a string type-name for the underlying vector. May appear in
	debug-output logs, etc. Does not need to exactly correspond to the type,
	but should be descriptive.
      */
      static const char* typeName();

      /** Set a specified scalar value throughout the vector.
       */
      static int setValues(T* vec, int firstLocalOffset,
			   double scalar, bool isSolnVector=false);

      /** Put values into the vector.
          If the 'sum_into' argument is true, then values are added to any
          values that already exist at the specified indices.
          If the 'sum_into' argument is false, then incoming values will
          overwrite any that may already exist at the specified indices.
          See general comments at the top of this struct regarding the
          'isSolnVector' argument.
      */
      static int putValuesIn(T* vec, int firstLocalOffset,
			     int numValues,
                             const int* indices,
                             const double* values,
                             bool sum_into,
                             bool isSolnVector=false,
                             int vectorIndex=0);

      /** Copy values for the specified indices out into the user-allocated
	  array 'values'.
      */
      static int copyOut(T* vec, int firstLocalOffset,
			 int numValues, const int* indices, double* values,
			 bool isSolnVector=false,
			 int vectorIndex=0);

      /** Get a pointer to the vector object's local coefficients array.
          Vector objects that can't support this can return NULL.
      */
      static double* getLocalCoefsPtr(T* vec,
                                      bool isSolnVector=false,
                                      int vectorIndex=0);

      /** Update vec = b*vec + a*x
       */
      static int update(T* vec,
			double a,
			const T* x,
			double b);

      /** Perform global communication or whatever operations may be
         necessary to complete the assembly of the vector. Most vector
         objects will do nothing here. Vectors such as the
         Epetra_FEVector object may do some operations here.
      */
      static int globalAssemble(T* vec);

    };//struct VectorTraits
}//namespace fei

#endif // _fei_VectorTraits_hpp_
