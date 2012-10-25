/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_FEVectorTraits_FED_hpp_
#define _snl_fei_FEVectorTraits_FED_hpp_

#include <fei_macros.hpp>

#include <snl_fei_FEVectorTraits.hpp>
#include <fei_FiniteElementData.hpp>

namespace snl_fei {

  /** Internal implementation vector traits. Define a "template" for accessing
      vector data.
  */
  template<>
  struct FEVectorTraits<FiniteElementData> {

    /** Return a string type-name for the vector. */
    static const char* typeName()
      { return("FiniteElementData"); }

    /** Reset (zero) the vector.
     */
    static int reset(FiniteElementData* vec)
      { return( vec->reset() ); }

    /** Sum an element-vector contribution into the FiniteElementData object */
    static int sumInElemVector(FiniteElementData* vec,
			       int elemBlockID,
			       int elemID,
			       int numNodes,
			       const int* nodeNumbers,
			       const int* dofPerNode,
             const int* dof_ids,
			       const double* coefs)
      {
	return( vec->setElemVector(elemBlockID, elemID, numNodes,
				   nodeNumbers, dofPerNode, dof_ids, coefs) );
      }

    /** Copy data out of the FiniteElementData object */
    static int copyOut(FiniteElementData* vec,
		       int nodeNumber,
		       int dofOffset,
		       double& value)
      {
	return( vec->getSolnEntry(nodeNumber, dofOffset, value) );
      }

  };//struct FEVectorTraits
}//namespace snl_fei

#endif // _snl_fei_FEVectorTraits_FED_hpp_
