/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef _snl_fei_FEVectorTraits_hpp_
#define _snl_fei_FEVectorTraits_hpp_

#include <fei_macros.hpp>

namespace snl_fei {

  /** Internal implementation vector traits. Define a "template" for accessing
      vector data.
      Provide function stubs for default type "T", which will catch the
      use of any vector type for which specialized traits have not been defined.
  */
  template<class T>
  struct FEVectorTraits {

    /** Return a string type-name for the vector. */
    static const char* typeName()
      { return("unsupported"); }

    /** Reset (zero) the vector.
     */
    static int reset(T* /*vec*/)
      { return(-1); }

    /** Sum an element-vector contribution into the vector object */
    static int sumInElemVector(T* /*vec*/,
			       int /*elemBlockID*/,
			       int /*elemID*/,
			       int /*numNodes*/,
			       const int* /*nodeNumbers*/,
			       const int* /*dofPerNode*/,
             const int* /*dof_ids*/,
			       const double* /*coefs*/)
      { return(-1); }

    /** Copy data out of the vector object */
    static int copyOut(T* /*vec*/,
		       int /*nodeNumber*/,
		       int /*dofOffset*/,
		       double& /*value*/)
      { return( -1 ); }

  };//struct FEVectorTraits
}//namespace snl_fei

#endif // _snl_fei_FEVectorTraits_hpp_
