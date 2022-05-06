/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef FaceSearchTolerance_hpp
#define FaceSearchTolerance_hpp

#include "stk_mesh/base/Types.hpp"

namespace stk {

namespace mesh { class BulkData; }
namespace mesh { class FieldBase; }
namespace mesh { struct Entity; }

namespace balance {

class FaceSearchTolerance
{
public:
  FaceSearchTolerance() {}
  virtual ~FaceSearchTolerance() {}

  virtual double compute(const stk::mesh::BulkData & mesh,
                         const stk::mesh::FieldBase & coordField,
                         const stk::mesh::Entity * faceNodes,
                         const unsigned numFaceNodes) const = 0;
};

}
}

#endif
