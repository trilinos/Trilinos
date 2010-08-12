/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_io_util_Skinning_hpp
#define stk_io_util_Skinning_hpp

namespace stk {
  namespace mesh {
    class BulkData;
    class Part;
  }

  namespace io {
    namespace util {

      void generate_sides(stk::mesh::BulkData & mesh,
			  stk::mesh::Part     &side_part,
			  const bool skin_only );
    } // namespace util
  } // namespace io
} // namespace stk
  
#endif 
