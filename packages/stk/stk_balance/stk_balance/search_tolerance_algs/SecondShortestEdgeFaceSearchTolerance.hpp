/*--------------------------------------------------------------------*/
/*    Copyright 2011 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SecondShortestEdgeFaceSearchTolerance_hpp
#define SecondShortestEdgeFaceSearchTolerance_hpp

#include "stk_balance/search_tolerance/FaceSearchTolerance.hpp"

namespace stk {
namespace balance {

class SecondShortestEdgeFaceSearchTolerance : public FaceSearchTolerance
{
public:
    SecondShortestEdgeFaceSearchTolerance()
      : FaceSearchTolerance() {}
    virtual ~SecondShortestEdgeFaceSearchTolerance() {}

    virtual double compute(const stk::mesh::BulkData & mesh,
                           const stk::mesh::FieldBase & coordField,
                           const stk::mesh::Entity * faceNodes,
                           const unsigned numFaceNodes) const override;

};

}
}

#endif
