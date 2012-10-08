/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_FEMatrixTraits_hpp_
#define _snl_fei_FEMatrixTraits_hpp_

#include <fei_macros.hpp>

namespace snl_fei {
  /** Internal implementation matrix traits. Define a "template" for accessing
      matrix data.
      Provide function stubs for default type "T", which will catch the
      use of any matrix type for which specialized traits have not been defined.
  */
  template<typename T>
  struct FEMatrixTraits {

    /** Return a string type-name for the underlying matrix */
    static const char* typeName()
      { return("unsupported"); }

    /** Reset (zero) the matrix.
     */
    static int reset(T* mat)
      { return(-1); }

    /** Sum an element-matrix contribution into the matrix object */
    static int sumInElemMatrix(T* mat,
			       int elemBlockID,
			       int elemID,
			       int numNodes,
			       const int* nodeNumbers,
			       const int* dofPerNode,
             const int* dof_ids,
			       const double *const * coefs)
      { return(-1); }

    /** Enforce Dirichlet boundary conditions on the matrix object */
    static int setDirichletBCs(T* mat,
			       int numBCs,
			       const int* nodeNumbers,
			       const int* dof_ids,
			       const double* values)
      { return(-1); }

  };//class FEMatrixTraits
}//namespace snl_fei

#endif // _snl_fei_FEMatrixTraits_hpp_
