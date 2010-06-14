/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef STK_MESH_UNITTEST_HEX_MESH_FIXTURE_HPP
#define STK_MESH_UNITTEST_HEX_MESH_FIXTURE_HPP

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/Types.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>

class HexFixture
{
public:
  HexFixture(stk::ParallelMachine pm);

  ~HexFixture();

  stk::mesh::MetaData& meta_data() { return m_meta_data; }
  stk::mesh::BulkData& bulk_data() { return m_bulk_data; }

  stk::mesh::Part* hex_part() const { return & m_hex_part; }
  stk::mesh::Part* skin_part() const { return & m_skin_part; }

private:
  void generate_mesh();

  stk::mesh::MetaData  m_meta_data;
  stk::mesh::BulkData  m_bulk_data;
  stk::mesh::Part    & m_hex_part;
  stk::mesh::Part    & m_skin_part;
};

#endif
