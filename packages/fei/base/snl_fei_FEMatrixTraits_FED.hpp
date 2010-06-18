/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_FEMatrixTraits_FED_hpp_
#define _snl_fei_FEMatrixTraits_FED_hpp_

#include <fei_macros.hpp>
#include <fei_FiniteElementData.hpp>
#include <snl_fei_FEMatrixTraits.hpp>

namespace snl_fei {

  /** specialization for FiniteElementData */
  template<>
  struct FEMatrixTraits<FiniteElementData> {

    /** Return a string type-name for the underlying matrix */
    static const char* typeName()
      { return("FiniteElementData"); }

    /** Reset (zero) the matrix.
     */
    static int reset(FiniteElementData* mat)
      { return( mat->reset() ); }

    /** sum-into operation for element-matrix data */
    static int sumInElemMatrix(FiniteElementData* mat,
			       int elemBlockID,
			       int elemID,
			       int numNodes,
			       const int* nodeNumbers,
			       const int* dofPerNode,
             const int* dof_ids,
			       const double *const * coefs)
      { return( mat->setElemMatrix(elemBlockID, elemID, numNodes,
				   nodeNumbers, dofPerNode, dof_ids, coefs) ); }

    /** specify dirichlet BCs */
    static int setDirichletBCs(FiniteElementData* mat,
			       int numBCs,
			       const int* nodeNumbers,
			       const int* dof_ids,
			       const double* values)
      { return( mat->setDirichletBCs(numBCs, nodeNumbers,
				     dof_ids, values) ); }

  };//struct FEMatrixTraits
}//namespace snl_fei

#endif // _snl_fei_FEMatrixTraits_FED_hpp_
