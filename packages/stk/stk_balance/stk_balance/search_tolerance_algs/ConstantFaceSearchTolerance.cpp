/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "ConstantFaceSearchTolerance.hpp"

namespace stk {
namespace balance {

double ConstantFaceSearchTolerance::compute(const stk::mesh::BulkData & /*mesh*/,
                                            const stk::mesh::FieldBase & /*coordField*/,
                                            const stk::mesh::Entity * /*faceNodes*/,
                                            const unsigned /*numFaceNodes*/) const
{
    return m_tolerance;
}

}
}
