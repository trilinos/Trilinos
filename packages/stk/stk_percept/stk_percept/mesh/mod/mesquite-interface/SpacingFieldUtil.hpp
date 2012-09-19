/*--------------------------------------------------------------------*/
/*    Copyright 2003 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef SpacingFieldUtil_hpp
#define SpacingFieldUtil_hpp

#include <stk_percept/Percept.hpp>
#if !defined(__IBMCPP__) && defined(STK_PERCEPT_HAS_MESQUITE)

#include <stk_percept/PerceptMesh.hpp>
#include "JacobianUtil.hpp"

namespace stk {
  namespace percept {

    class SpacingFieldUtil 
    {
    public:
      /// either average or take max of spacing to the nodes
      enum SpacingType { SPACING_AVE, SPACING_MAX };

      SpacingFieldUtil(PerceptMesh& eMesh, SpacingType type=SPACING_AVE) : m_eMesh(eMesh), m_type(type) {}


      void compute_spacing_field();

    private:
      PerceptMesh& m_eMesh;
      SpacingType m_type;
      
    };
  }
}

#endif
#endif
